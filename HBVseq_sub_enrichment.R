# Ashley Acevedo
# 2018-07-16
# Last update: 2019-02-06
# Step 1: Read reference sequence
# Step 2: Extract codons from BAM files, generate count tables and calculate statistics for each substitution
# Step 3: Combine p.values from multiple experiments (if multiple eperimental sets are provided)
# Step 4: Write mutant data to a file
# Step 5: Plot 1: ratio of frequencies versus sum of counts (confidence) for each set
# Step 6: Plot 2: -log10 p.value versus residue position

library(Rsamtools)
library(seqinr)
library(metap)
library(data.table)

################################################################################
################################### Inputs: ####################################
################################################################################

# Run date
run.date = "2019-02-06"

# Define the sets to be compared. Add as many sets as desired as long as the set
# has a "name" (below). Set "a" must be used.

# Set names: 
sets = c("a", "b", "c")

# A
set.a.BAM.PATH = "bam_files/"
set.a.BAM.FILEs = c("Set_A_drug.bam", "Set_A_no_drug.bam")

# B
set.b.BAM.PATH = "bam_files/"
set.b.BAM.FILEs = c("Set_B_drug.bam", "Set_B_no_drug.bam")

# C
set.c.BAM.PATH = "bam_files/"
set.c.BAM.FILEs = c("Set_C_drug.bam", "Set_C_no_drug.bam")

# Reference sequence
REFERENCE.PATH = "sample_reference/"
REFERENCE = "HBV_amplicon.fasta"

# For multiple sets, choose a method for combining p values
#  options: "fisher", "z.transform" or "weighted.z.transform"
#  For a single set, fisher's exact test p values are reported
method = "weighted.z.transform"

# Define output path and filenames/prefix
output.PATH       = "outputs/"
top.hits.filename = "Top_substitutions.txt"
all.data.filename = "All_substitutions.txt"
plot.1.prefix     = "log10_counts_v_ratio_set_"
plot.2.name       = "log10_p_v_substitution"

# Enter TRUE to produce output, FALSE otherwise
output.top.hits = TRUE
output.all.data = TRUE
output.plot.1   = TRUE
output.plot.2   = TRUE

# Number of hits to output
HITS = 20

# Enter the range of amino acid positions to output (use to exclude regions 
# with poor expected data quality). For HBV, 1:356 is the full RT sequence. Any 
# vector of numbers in any order will work in this argument.
output.range = 1:356

# Save workspace for later data manipulation
save.workspace = TRUE
workspace.name = "_HBVseq_run"

################################################################################
########################### User-defined Functions: ############################
################################################################################

# Function to: Break reads into codons.
#			   Ajusts frame based on read start position.
#			   Ignores reads with indels and low quality bases.
#			   Renumbers codons based on amino acid sequence.

SplitCodons = function(sequence, quality, position, cigar){
	
	cigar = strsplit(cigar, split = "", fixed = TRUE)[[1]]
	
	if ("I" %in% cigar | "D" %in% cigar | "N" %in% cigar | anyNA(cigar)){
		list("codons"    = NA,
			 "positions" = NA)
	} else {
		
		# Find/adjust frame of each sequence
		seq.mod  = position %% 3
		sequence = as.character(sequence)
		quality  = as.character(quality)
		
		if (seq.mod == 2){
			sequence = substring(sequence, 3)
			quality  = substring(quality , 3)
			position = position + 2
		} else if (seq.mod == 0){
			sequence = substring(sequence, 2)
			quality  = substring(quality , 2)
			position = position + 1			
		}
		
		# Break by codons
		codons.seq = strsplit(sequence, "(?<=.{3})", perl = TRUE)[[1]]
		codons.q   = strsplit(quality , "(?<=.{3})", perl = TRUE)[[1]]
		codons.pos = seq.int(from = position, by = 3, length.out = length(codons.seq))		

		# Ignore codons with < 3 or low quality bases
		codons    = c()
		positions = c()
		for (j in 1:length(codons.seq)){
			if (nchar(codons.q[j]) == 3){
				if(min(utf8ToInt(codons.q[j]) - 33) >= 13){
					codons    = c(codons,    codons.seq[j])
					positions = c(positions, codons.pos[j])
				}
			}
		}
		
		if (length(codons) > 0){
			
			# Renumber codons to make sequential (and with respect to amino acid sequence)
			positions = (positions - 1)/3 + 1
			
			list("codons"    = codons,
		     	 "positions" = positions)
		} else {
			list("codons"    = NA,
			 	 "positions" = NA)
		}
	}
}

# Function to: Generate amino acid count table.

MakeCountTable = function(codon.list){
	
	
	# Initialize count.table
	count.table = matrix(0, ncol = 21, nrow = no.codons)
	count.table = as.data.frame(count.table)
	colnames(count.table) = c("A", "C", "D", "E", "F", "G", "H",
							  "I", "K", "L", "M", "N", "P", "Q",
							  "R", "S", "T", "V", "W", "Y", "*")
	
	# Convert list of codons and positions into a data frame
	codon.df = lapply(codon.list, data.frame, stringsAsFactors = FALSE)
	codon.df = rbindlist(codon.df)
	codon.df = codon.df[complete.cases(codon.df), ]
	
	# Translate codons and populate count.table with amino acid counts
	for (i in 1:max(codon.df$positions)){
		codon.table = table(codon.df$codons[codon.df$positions == i])
		amino.acids = sapply(rownames(codon.table), TranslateCodons)
		codon.table = as.data.frame(codon.table) # Column names default to Var1 and Freq (for codon and counts)
		codon.table$amino.acids = amino.acids
		for(j in 1:nrow(codon.table)){
			count.table[[codon.table$amino.acids[j]]][i] = count.table[[codon.table$amino.acids[j]]][i] + codon.table$Freq[j]
		}
	}
	
	# Add reference amino acid sequence and position to count.table
	count.table$reference = reference.trans
	count.table$position  = 1:length(reference.trans) - 13
	
	return(count.table)				  
}

# Function to: Translate codons.

TranslateCodons = function(codon){
	codon = strsplit(codon, split = "", fixed = TRUE)[[1]]
	amino.acid = translate(codon, frame = 0)
	return(amino.acid)
}

# Function to: Perform Fisher's exact test

FishersExactTest = function(mutant.cond.A, total.cond.A, mutant.cond.B, total.cond.B){
	
	contingency.table = matrix(c(mutant.cond.A, mutant.cond.B, total.cond.A - mutant.cond.A, total.cond.B - mutant.cond.B), nrow = 2)
	fisher.result = fisher.test(contingency.table, alternative = "greater")
	return(fisher.result$p.value)
	
}

# Function to: Run weighted z-transform test

WeightedZTransformTest = function(p.and.n){
	
	no.tests = length(p.and.n)/2
	WEIGHTS = sqrt(p.and.n[no.tests + 1: no.tests])
	return(sumz(p.and.n[1:no.tests], weights = WEIGHTS)$p)
	
}

################################################################################
#################################### Main: #####################################
################################################################################

# Step 1: Read reference sequence.

reference       = read.fasta(file = paste(REFERENCE.PATH, REFERENCE, sep = ""), seqtype = "DNA", forceDNAtolower = FALSE)[[1]]
reference.trans = seqinr::translate(reference)
no.codons       = length(reference.trans)

# Step 2: For each set: Extract codons from BAM files,
#						Generate count tables,
#						Determine count frequencies for each substitution,
#						Perform fisher's exact test on counts for each mutation,
#						Compile a dataframe with outputs for further tests/plots

for (i in sets){
	
	# From inputs
	BAM.PATH = get(paste("set.", i, ".BAM.PATH", sep = ""))
	BAM.FILEs = get(paste("set.", i, ".BAM.FILEs", sep = ""))
	
	# Read .BAM files and extract codons
	print("reading file 1")
	bam.data = scanBam(paste(BAM.PATH, BAM.FILEs[1], sep = ""))
	drug.treated.codon.list = mapply(function(w, x, y, z) SplitCodons(w, x, y, z),
									 bam.data[[1]]$seq,
									 bam.data[[1]]$qual,
									 bam.data[[1]]$pos,
									 bam.data[[1]]$cigar,
									 SIMPLIFY = FALSE)
									 
	print("reading file 2")
	bam.data = scanBam(paste(BAM.PATH, BAM.FILEs[2], sep = ""))								 
	drug.untreated.codon.list = mapply(function(w, x, y, z) SplitCodons(w, x, y, z),
									 bam.data[[1]]$seq,
									 bam.data[[1]]$qual,
									 bam.data[[1]]$pos,
									 bam.data[[1]]$cigar,
									 SIMPLIFY = FALSE)	

	# Convert lists of codons into amino acid count tables
	print("making count table 1")
	drug.treated.count.table   = MakeCountTable(drug.treated.codon.list)
	print("making count table 2")
	drug.untreated.count.table = MakeCountTable(drug.untreated.codon.list)

	# Initialize a dataframe to populate below
	temp.df = data.frame(ref.amino.acid 				  = character(),
						 mutant.amino.acid 			  = character(),
						 position 					  = numeric(),
						 ratio.of.freq 				  = numeric(),
						 pseudo.ratio.of.freq		  = numeric(),
						 sum.of.counts 				  = numeric(),
						 n.counts 					  = numeric(),
						 p.values 					  = numeric(),
						 drug.treated.mutant.counts   = numeric(),
						 drug.untreated.mutant.counts = numeric(),
						 drug.treated.total.counts    = numeric(),
						 drug.untreated.total.counts  = numeric())

	# For each substitution type: Calculate count frequency in each condition
	#						  	  Calculate the ratio of frequencies
	#						  	  Calculate a confidence measure: sum.of.counts
	#						  	  Perform Fisher's exact test on counts
	#						      Append data to temp.result.df
	
	amino.acids = c("A", "C", "D", "E", "F", "G", "H",
					"I", "K", "L", "M", "N", "P", "Q",
					"R", "S", "T", "V", "W", "Y", "*")
	
	print("analyzing results")				
	for (a in amino.acids){
		
		# Subset count table by reference amino acid
		drug.treated.count.table.subset = drug.treated.count.table[which(drug.treated.count.table$reference == a), ]
		drug.untreated.count.table.subset = drug.untreated.count.table[which(drug.untreated.count.table$reference == a), ]

		# For each substitution type: Perform calculations noted above
		for (s in amino.acids){
			
			if (s != a){
				
				# Calculate count frequency
				drug.treated.mutant.counts  = drug.treated.count.table.subset[[s]]
				drug.untreated.mutant.counts = drug.untreated.count.table.subset[[s]]
				drug.treated.total.counts    = rowSums(drug.treated.count.table.subset[ , 1:21])
				drug.untreated.total.counts  = rowSums(drug.untreated.count.table.subset[ , 1:21])
				drug.treated.freq   = drug.treated.mutant.counts/drug.treated.total.counts
				drug.untreated.freq = drug.untreated.mutant.counts/drug.untreated.total.counts
				pseudo.drug.untreated.freq = 1/rowSums(drug.untreated.count.table.subset[ , 1:21])
			
				# Calulate ratio of frequencies
				ratio.of.freq = drug.treated.freq/drug.untreated.freq
			
				# Calculate "confidence": sum of counts
				sum.of.counts = drug.treated.count.table.subset[[s]] + drug.untreated.count.table.subset[[s]]
			
				# Calculate pseudo frequencies (for no count untreated positions)
				pseudo.ratio.of.freq = drug.treated.freq/pseudo.drug.untreated.freq
				
				# Fisher's exact test for count data
				p.values = mapply(function(a, b, c, d) FishersExactTest(a, b, c, d), 
							  		   drug.treated.mutant.counts,   
							  		   drug.treated.total.counts,
									   drug.untreated.mutant.counts, 
									   drug.untreated.total.counts,
									   SIMPLIFY = TRUE)

				# Combine data into dataframe
				no.mutants  		  = length(ratio.of.freq)
				ref.amino.acid    = rep(a, no.mutants)
				mutant.amino.acid = rep(s, no.mutants)
				position    		  = drug.treated.count.table.subset$position
				n.counts    		  = drug.treated.total.counts + drug.untreated.total.counts
				mutant.df   		  = data.frame(ref.amino.acid, 
									 	 	   mutant.amino.acid, 
									 	 	   position, 
									 	 	   ratio.of.freq, 
									 	 	   pseudo.ratio.of.freq, 
									 	 	   sum.of.counts,
									 	 	   n.counts,
									 	 	   p.values,
									 	 	   drug.treated.mutant.counts,
									 	 	   drug.untreated.mutant.counts,
									 	 	   drug.treated.total.counts,
									 	 	   drug.untreated.total.counts)
			
				# Append ratio.df to the dataframe for the entire set
				temp.df = rbind(temp.df, mutant.df)
				
			}
		}
	}
	
	assign(paste("set.", i, ".df", sep = ""), temp.df)
}

# Step 3: Combine p.values from multiple, independent experiments (if multiple sets are provided)
# For reference: Whitlock, MC. Combining probability from independent tests:
#				 	the weighted Z-method is superior to Fisher's approach.
#					(2005) J Evol Biol. 18: 1368-1373.

if (length(sets) > 1){

	# Empty vectors to populate with p.values and n from each experiment
	p.value.df = c()
	n.df = c()

	# Extract p.values and n from each set and append to above vectors
	for (i in sets){
		
		temp.df = get(paste("set.", i, ".df", sep = ""))
		p.value.df = cbind(p.value.df, temp.df$p.value)
		n.df = cbind(n.df, temp.df$n.counts)
	
	}

	# Alter p.value.df to eliminate p = 1, which is incompatible with z-transform tests
	alt.p.value.df = p.value.df
	alt.p.value.df[alt.p.value.df > 0.9999999] = 0.9999999

	# Combine p.values using 3 methods
	fishers.combined.p.values = apply(p.value.df, MARGIN = 1, function(x) sumlog(x)$p)
	z.transform.p.values = apply(alt.p.value.df, MARGIN = 1, function(x) sumz(x)$p)
	weighted.z.transform.p.values = apply(cbind(alt.p.value.df, n.df), MARGIN = 1, function(x) WeightedZTransformTest(x))

	# Generate table of combined p.values for sorting top hits
	combined.p.values.df = data.frame(set.a.df$ref.amino.acid, set.a.df$mutant.amino.acid, set.a.df$position, fishers.combined.p.values, z.transform.p.values, weighted.z.transform.p.values)
	colnames(combined.p.values.df) = c("ref.amino.acid", "substitution", "position","fisher","z.transform", "weighted.z.transform")

	# Add ratio of counts, sum of counts and raw counts to combined.p.values.df
	for (i in sets){
	
		temp.df = get(paste("set.", i, ".df", sep = ""))
		combined.p.values.df[[paste("set.", i, ".ratio.of.freq", sep = "")]] 		   =  temp.df$ratio.of.freq
		combined.p.values.df[[paste("set.", i, ".sum.of.counts", sep = "")]] 		   =  temp.df$sum.of.counts
		combined.p.values.df[[paste("set.", i, ".treated.mutant.counts", sep = "")]]   =  temp.df$drug.treated.mutant.counts
		combined.p.values.df[[paste("set.", i, ".untreated.mutant.counts", sep = "")]] =  temp.df$drug.untreated.mutant.counts
		combined.p.values.df[[paste("set.", i, ".treated.total.counts", sep = "")]]    =  temp.df$drug.treated.total.counts
		combined.p.values.df[[paste("set.", i, ".untreated.total.counts", sep = "")]]  =  temp.df$drug.untreated.total.counts
	
	}
		
} else if(length(sets) == 1){
	
	# Generate table of p.value data for sorting top hits
	temp.df = get(paste("set.", sets[1], ".df", sep = ""))
	p.values.df = data.frame(temp.df$ref.amino.acid, 
							 temp.df$mutant.amino.acid, 
							 temp.df$position, 
							 temp.df$p.values, 
							 temp.df$ratio.of.freq, 
							 temp.df$sum.of.counts,
							 temp.df$drug.treated.mutant.counts,
							 temp.df$drug.untreated.mutant.counts,
							 temp.df$drug.treated.total.counts,
							 temp.df$drug.untreated.total.counts)
	colnames(p.values.df) = c("ref.amino.acid", 
							  "substitution", 
							  "position", 
							  "fisher", 
							  paste("set.", sets[1], ".ratio.of.freq", sep = ""), 
							  paste("set.", sets[1], ".sum.of.counts", sep = ""),
							  paste("set.", sets[1], ".treated.mutant.counts", sep = ""),
							  paste("set.", sets[1], ".untreated.mutant.counts", sep = ""),
							  paste("set.", sets[1], ".treated.total.counts", sep = ""),
							  paste("set.", sets[1], ".untreated.total.counts", sep = ""))

}

################################################################################
################################## Outputs: ####################################
################################################################################

# Step 4: Write tables of mutations by p.value

# Sort data by p.value
if (length(sets) > 1){
	sort.by = method
	sorted.df = combined.p.values.df[order(combined.p.values.df[[sort.by]]), ]
} else if (length(sets) == 1){
	sort.by = "fisher"
	sorted.df = p.values.df[order(p.values.df[[sort.by]]), ]
}

# Filter out positions outside of the range set by output.range
sorted.df = sorted.df[sorted.df$position %in% output.range, ]

# Build a subheader for ratio of frequency and sum of counts dependent on the number of sets input
header.string.of.sets = ""

for (i in sets){
	
	header.string.of.sets = paste(header.string.of.sets, "ratio.of.freq.set_",   		   i, 
														 "\tsum.of.substitution.counts.set_", 		   i, 
														 "\ttreated.substitution.counts.set_",   i,
														 "\ttreated.total.counts.set_",    i,
														 "\tuntreated.substitution.counts.set_", i,
														 "\tuntreated.total.counts.set_",  i,
														 "\t", sep = "")
	
}

if (output.all.data){

	# Write file header
	cat(paste("Reference", "Position", "Substitution", "p.value", header.string.of.sets, sep = "\t"), 
    		file = paste(output.PATH, all.data.filename, sep = ""), 
    		sep = "\n", 
    		append = FALSE)

	# Write data
	for (i in 1:nrow(sorted.df)){
	
		set.data = ""
		for (j in sets){
		
			set.data = paste(set.data, sorted.df[[paste("set.", j, ".ratio.of.freq", sep = "")]][i], "\t", 
									   sorted.df[[paste("set.", j, ".sum.of.substitution.counts", sep = "")]][i], "\t", 
									   sorted.df[[paste("set.", j, ".treated.substitution.counts", sep = "")]][i], "\t",
									   sorted.df[[paste("set.", j, ".treated.total.counts", sep = "")]][i], "\t",
									   sorted.df[[paste("set.", j, ".untreated.substitution.counts", sep = "")]][i], "\t",
									   sorted.df[[paste("set.", j, ".untreated.total.counts", sep = "")]][i], "\t",
									   sep = "")
		}

		cat(paste(sorted.df$ref.amino.acid[i], "\t",
			  	  sorted.df$position[i], "\t",
			  	  sorted.df$substitution[i], "\t",
			  	  signif(sorted.df[[sort.by]][i], 3), "\t",
			  	  set.data, sep = ""),
			  	  file = paste(output.PATH, all.data.filename, sep = ""),
			   	  sep = "\n",
			  	  append = TRUE )
	}
}

if (output.top.hits){

	# Write header
	cat(paste("Reference", "Position", "Substitution", "p.value", header.string.of.sets, sep = "\t"), 
    		file = paste(output.PATH, top.hits.filename, sep = ""), 
    		sep = "\n", 
    		append = FALSE)
    
	# Write data
	for (i in 1:nrow(sorted.df)){
	
		set.data = ""
		for (j in sets){
		
			set.data = paste(set.data, sorted.df[[paste("set.", j, ".ratio.of.freq", sep = "")]][i], "\t", 
									   sorted.df[[paste("set.", j, ".sum.of.substitution.counts", sep = "")]][i], "\t", 
									   sorted.df[[paste("set.", j, ".treated.substitution.counts", sep = "")]][i], "\t",
									   sorted.df[[paste("set.", j, ".treated.total.counts", sep = "")]][i], "\t",
									   sorted.df[[paste("set.", j, ".untreated.substitution.counts", sep = "")]][i], "\t",
									   sorted.df[[paste("set.", j, ".untreated.total.counts", sep = "")]][i], "\t",
									   sep = "")
		}
			  
		if (i <= HITS){
		
			cat(paste(sorted.df$ref.amino.acid[i], "\t",
			  	  	  sorted.df$position[i], "\t",
			  	  	  sorted.df$substitution[i], "\t",
			  	  	  signif(sorted.df[[sort.by]][i], 3), "\t",
			  	  	  set.data, sep = ""),
			  	  	  file = paste(output.PATH, top.hits.filename, sep = ""),
			   	  	  sep = "\n",
			  	  	  append = TRUE  )
		
		} else { break }
	}
}

# Step 5: Plot ratio.of.freq vs sum.of.counts for each set
if (output.plot.1){
	for (i in sets){
		
		temp.df = get(paste("set.", i, ".df", sep = ""))
		pseudo.df = temp.df[which(temp.df$ratio.of.freq == Inf), ]
		significant.df = temp.df[which(temp.df$p.values <= 0.001), ]
		pseudo.significant.df = pseudo.df[which(pseudo.df$p.values <= 0.001), ]
	
		pdf(file = paste(output.PATH, plot.1.prefix, i, ".pdf", sep = ""), width = 5, height = 5)
		plot(log10(temp.df$sum.of.counts), temp.df$ratio.of.freq, 
			 ylim = c(0, 12), 
			 xlim = c(0, 3), 
			 xlab = "log10(Sum of substitution counts)", 
			 ylab = "Ratio of substitution frequencies", 
			 pch  = 16, 
			 col  = "grey", 
			 axes  = FALSE)
		points(log10(pseudo.df$sum.of.counts), pseudo.df$pseudo.ratio.of.freq,
			 pch = 16,
			 col = "yellow")
		points(log10(significant.df$sum.of.counts), significant.df$ratio.of.freq,
			 pch = 16,
			 col = "black")
		points(log10(pseudo.significant.df$sum.of.counts), pseudo.significant.df$pseudo.ratio.of.freq,
			 pch = 16,
			 col = "gold4")													  	
		axis(1, at = seq(0, 3, 0.5), labels = seq(0, 3, 0.5))
		axis(2, at = seq(0, 12, 2), labels = seq(0, 12, 2), las = 1)
		dev.off()
	
	}
}	

# Step 6: Plot -log10 p.value versus residue position

# Remove substitutions not observed in any sample
substitution.counts = rep(0, nrow(sorted.df))
for (i in sets){	
	substitution.counts = substitution.counts +  sorted.df[[paste("set.", i, ".sum.of.counts", sep = "")]]
}
sorted.df.culled = sorted.df[substitution.counts > 0, ]

# Save the workspace
if (save.workspace == TRUE){save.image(paste(output.PATH, run.date, workspace.name, ".RData", sep = ""))}

if (output.plot.2){

	pdf(file = paste(output.PATH, plot.2.name, ".pdf", sep = ""), height = 4, width = 6)
	plot(sorted.df.culled$position, -log10(sorted.df.culled[[sort.by]]), 
		 pch = 20,
		 xlab = "Residue number",
		 ylab = expression("-log"[10]*"(p-value)"),
		 ylim = c(0, 125),
		 col = "grey",
		 axes = FALSE)
	axis(1, at = seq(0, 350, 50), labels = seq(0, 350, 50))
	axis(2, at = seq(0, 125, 25), labels = seq(0, 125, 25), las = 1)

	# Add line to note 0.05 bonferroni correction threshold
	abline(h = -log10(0.05/nrow(sorted.df.culled)), lwd = 2, lty = 2, col = "grey")

	# Positive hits
	hit.positions    = c(204, 204)
	hit.substitution = c("V", "I")
	hit.name         = c("M204V", "M204I")
	for (i in 1:length(hit.positions)){
		points(sorted.df.culled$position[sorted.df.culled$position == hit.positions[i] & sorted.df.culled$substitution == hit.substitution[i]], 
			   -log10(sorted.df.culled[[sort.by]][sorted.df.culled$position == hit.positions[i] & sorted.df.culled$substitution == hit.substitution[i]]), 
			   col = "dodgerblue", pch = 20)

		text(sorted.df.culled$position[sorted.df.culled$position == hit.positions[i] & sorted.df.culled$substitution == hit.substitution[i]], 
			 -log10(sorted.df.culled[[sort.by]][sorted.df.culled$position == hit.positions[i] & sorted.df.culled$substitution == hit.substitution[i]]), 
			 col = "dodgerblue", labels = hit.name[i], pos = 4, cex = 0.5)

	}

	# Other substitutions to highlight
	oth.positions    = c(80, 180, 181, 181, 181, 194, 173, 236)
	oth.substitution = c("I", "M", "T", "V", "S", "T", "L", "T")
	for (i in 1:length(oth.positions)){
		points(sorted.df.culled$position[sorted.df.culled$position == oth.positions[i] & sorted.df.culled$substitution == oth.substitution[i]], 
	  		   -log10(sorted.df.culled[[sort.by]][sorted.df.culled$position == oth.positions[i] & sorted.df.culled$substitution == oth.substitution[i]]), 
	  		   col = "dimgrey", pch = 20)
	}

	dev.off()
}