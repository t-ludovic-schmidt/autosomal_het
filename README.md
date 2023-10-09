# Autosomal Heterozygosity
For deriving heterozygosity estimates from genomic data without filtering for polymorphism.

This code was developed for the paper: 
### Schmidt TL, Thia JA, and Hoffmann AA. (2024) How Can Genomics Help or Hinder Wildlife Conservation? Annual Review of Animal Biosciences. 12 (1). <br>https://doi.org/10.1146/annurev-animal-021022-051810 
<br> <br>
If you use this code, please cite the above paper. <br> <br>
This code builds on previous work (Schmidt et al. 2021 Meth Ecol Evol; https://doi.org/10.1111/2041-210X.13659). Optimal heterozygosity workflows will retain all individuals at the genotyping stage, but then filter each individual separately for depth and missing data, and will retain all monomorphic, biallelic, and polyallelic sites that pass these filters. This is one such workflow.

This code takes as inputs a set of .bam files (one per individual) aligned to the same reference, and outputs a set of text files (one per individual) containing genotype calls at polymorphic and monomorphic sites. The reference genome can be of any quality, even highly fragmented genomes of 10,000s - 100,000s of contigs are appropriate, though additional parameters in GATK may be needed to decrease run time. <br> <br>
The early filtering and HaplotypeCaller steps can be run on subsets of samples as needed, but the GenotypeGVCFs step should be run on all individuals together. If dataset is too large, analyse subsets of intervals, build each subset into its own database, genotype and filter each subset individually, then merge filtered VCFs into a final VCF. <br> <br>
This code has been developed using GATK v4.2.6.1, Picard v2.27.4, samtools v1.16, and bcftools v1.16 
