HBVseq_sub_enrichment
=======================
Description
---------------------------------------------------------------------------
HBVseq_sub_enrichment is an application designed to identify amino acid 
substitutions enriched in virus populations following drug treatment of 
infected cells. As input, the application requires:

1. Paired sets of high depth sequence aligned to a reference in BAM format 

	Each set should consist of one "treated" and one "untreated"
	sample. All sets must be aligned to the same reference sequence.

2. A reference sequence in FASTA format
	
	The reference sequence should be the same as that used to generate
	the BAM formatted alignments. The reference should be in-frame.

The following outputs are generated:

1. Top_substitutions.txt

	Contains a tab-separated table of the top substitutions sorted by
	p-value (low to high). Each substitution is numbered based on its
	position in the reference sequence, with the first codon in the
	reference being position 1. For a single set, the p-value is from
	a fisher's exact test. For multiple sets, the p-value is combined
	for all the sets using either fisher's method, a z-transform or a 
	weighted z-transform. Weighted z-transform is the default but can 
	be changed in the INPUTS section (line 48). For each set, the table 
	contains counts of total codons and substitutions at each site in 
	both "treated" and "untreated" conditions, as well as the sum of 
	substitution counts in both conditions and the ratio of 
	substitution frequency ("treated"/"untreated").

2. All_substitutions.txt

	Contains a tab-separated table of all substitutions sorted by
	p-value (low to high) as described above.

3. log10_counts_v_ratio_set_*.pdf
	
	For each paired sample set, a scatter plot of the log10(sum of 
	substitution counts) versus ratio of substitution frequency for 
	each substitution is generated. Data points in black are 
	substitutions with p-values <= 0.001. Data points in yellow are
	based on a pseudo count. For cases in which a substitution is 
	observed in the "treated" condition but not in the "untreated"
	condition, the ratio of substitution frequencies is infinity. In
	this case, a single substitution count is added to the "untreated"
	sample, a pseudo count, to support plotting. Gold points are
	substitutions with a p-value <= 0.001 that required a pseudo count 
	for plotting. 

4. log10_p_v_substitution.pdf

	A scatter plot of log10(combined p-value) versus substitution 
	position is generated. Horizontal dashed line corresponds to 0.05
	bonferroni correction threshold.

System requirements
---------------------------------------------------------------------------
The following packages are prerequisites for running HBVseq_sub_enrichment:

1. R (version 3.3.3)
2. R packages: Rsamtools, seqinr, metap, data.table

Usage
---------------------------------------------------------------------------
Prior to running, modify the parameters set in the INPUTS sections on lines
20-73. The appearance of plots should be modified to accommodate the range
of values obtained for each run. Plot 1 can be modified on lines 528-546.
Plot 2 can be modified on lines 567-575. Substitutions of interest can be
highlighted by adding to lists on lines 581-583 and 596-597.

On command line:

	Rscript HBVseq_sub_enrichment.R

In interactive console:

	source("HBVseq_sub_enrichment.R")

Test set
---------------------------------------------------------------------------
Sample data is provided (by Yingpu Yu, Ph.D.) in the bam_files directory.
This sample data includes three paired sets of treated and untreated virus
sequence aligned to the reference included in the sample_reference
directory.
