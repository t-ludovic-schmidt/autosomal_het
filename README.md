# Autosomal Heterozygosity
Code for deriving heterozygosity estimates from genomic data without filtering for polymorphism.

This code was developed for the paper: 
### Schmidt TL, Thia JA, and Hoffmann AA. (2023) How Can Genomics Help or Hinder Wildlife Conservation? Annual Review of Animal Biosciences 
If you use this code, please cite this paper.

This code takes as inputs a set of .bam files aligned to the same reference, and outputs a set of text files containing genotype calls at polymorphic and monomorphic sites. The reference genome can be of any contiguity, even highly fragmented genomes of 10,000s - 100,000s of contigs are appropriate, though additional parameters in GATK may be needed to decrease run time. <br> <br>
The early filtering and HaplotypeCaller steps can be run on subsets of samples as needed, but the GenotypeGVCFs step should be run on all individuals together. If dataset is too large, analyse subsets of intervals, build each subset into its own database, genotype and filter each subset individually, then merge filtered VCFs into a final VCF. <br> <br>
This code has been developed using GATK v4.2.6.1, Picard v2.27.4, samtools v1.16, and bcftools v.1.16 
