#!/bin/bash

#### README
#### This code was developed for the paper: "Schmidt TL, Thia JA, and Hoffmann AA. (2023) How Can Genomics Help or Hinder Wildlife Conservation? Annual Review of Animal Biosciences". If you use this code, please cite this paper.
#### This code takes as inputs a set of .bam files aligned to the same reference, and outputs a set of text files containing genotype calls at polymorphic and monomorphic sites.
#### The genome can be of any contiguity, even highly fragmented genomes of 10,000s - 100,000s of contigs are appropriate, though additional parameters in GATK may be needed to decrease run time.
#### The early filtering and HaplotypeCaller steps can be run on subsets of samples as needed, but the GenotypeGVCFs step should be run on all individuals together. If dataset is too large, analyse subsets of intervals, build each subset into its own database, genotype and filter each subset individually, then merge filtered VCFs into a final VCF.
#### This code has been developed using GATK v4.2.6.1, Picard v2.27.4, samtools v1.16, and bcftools v.1.16
#### To run this code, first set all paths as below


#### Define variables
#### See GATK documentation for format of GenomicsDB sample map and intervals files
samples="Sample1 Sample2 Sample3 Sample4 Sample5"
gatk="PATH/TO/GATK"
picard="PATH/TO/PICARD"
ref="PATH/TO/REFERENCE.fasta|.fa"
bam_dir="PATH/TO/INPUT_BAM_DIR"
workspace="PATH/TO/OUTPUT_FILES_DIR"
db_path="PATH/TO/WRITE_DATABASE"
db_samplemap="PATH/TO/MAP.sample_map"
intervals="PATH/TO/CHR.intervals"



#### Create reference dict and fai files
$gatk CreateSequenceDictionary -R $ref
samtools faidx $ref


#### Filter bam files
#### If using WGS or similar data, run Picard's MarkDuplicates step. DO NOT run MarkDuplicates with ddRAD or similar data
#### If bam files do not already have read groups, run Picard's AddOrReplaceReadGroups to assign
for K in $samples;
do
	java -jar $picard AddOrReplaceReadGroups \
		I=$bam_dir/${K}.bam \
		O=$bam_dir/${K}_rg.bam \
		SORT_ORDER=coordinate \
		RGID=1 \
		RGLB=lib1 \
		RGPL=illumina \
		RGPU=unit1  \
		RGSM=${K} \
		CREATE_INDEX=True
	#java -jar $picard MarkDuplicates \
	#	I=$bam_dir/${K}_rg.bam \
	#	O=$bam_dir/${K}_md.bam \
	#	M=$bam_dir/${K}_dup.txt
	java -jar $picard ValidateSamFile \
		I=$bam_dir/${K}_rg.bam \
		MODE=SUMMARY
	samtools index $bam_dir/${K}_rg.bam
	samtools view -F 4 -F 256 -f 0x02 -h $bam_dir/${K}_rg.bam | \
		samtools view -b -T $ref - > $bam_dir/${K}_filt.bam
	samtools index $bam_dir/${K}_filt.bam
done


#### Produce gVCFs for each sample 
for K in $samples;
do
	$gatk --java-options "-Xmx16g" HaplotypeCaller  \
		-R $ref  \
		-I $bam_dir/${K}_filt.bam \
		--native-pair-hmm-threads 4  \
		-output-mode EMIT_ALL_CONFIDENT_SITES  \
		-ERC GVCF \
		--sample-ploidy 2  \
		-O $workspace/${K}_g.vcf.gz
	$gatk --java-options "-Xmx16g" IndexFeatureFile \
		-I $workspace/${K}_g.vcf.gz
done


#### Build database. 
#### If an existing database has the same name, delete or rename
$gatk --java-options "-Xmx16g" GenomicsDBImport \
	--sample-name-map $db_samplemap \
	--genomicsdb-workspace-path $db_path \
	--merge-contigs-into-num-partitions 3 \
	--L $intervals


#### Genotype all samples together. 
#### If memory limits are reached, reduce number or size of intervals in each run but retain all individuals. 
#### Large chromosomes may need to be analysed as sections. 
#### See concat function after hard filtering for merging filtered VCFs into one whole.
$gatk --java-options "-Xmx16g" GenotypeGVCFs \
	-R $ref  \
	-V gendb://$db_path  \
	--add-output-vcf-command-line true  \
	--include-non-variant-sites true  \
	--sample-ploidy 2  \
	--L $intervals \
	-O $workspace/allgeno.vcf.gz


#### Remove indels and apply hard filtering 
$gatk --java-options "-Xmx16g" SelectVariants \
	-V $workspace/allgeno.vcf.gz \
	--select-type-to-exclude INDEL \
	--verbosity ERROR \
	-O $workspace/noindel.vcf.gz
	
$gatk --java-options "-Xmx16g" VariantFiltration \
	-V $workspace/noindel.vcf.gz \
	--verbosity ERROR \
	-filter "QD < 2.0" --filter-name "QD2" \
	-filter "QUAL < 30.0" --filter-name "QUAL30" \
	-filter "SOR > 3.0" --filter-name "SOR3" \
	-filter "FS > 60.0" --filter-name "FS60" \
	-filter "MQ < 40.0" --filter-name "MQ40" \
	-O $workspace/filt.vcf.gz
	
$gatk --java-options "-Xmx16g" SelectVariants \
	-V $workspace/filt.vcf.gz \
	--verbosity ERROR \
	--set-filtered-gt-to-nocall TRUE \
	-O $workspace/final.vcf.gz


#### OPTIONAL: If you have run genotyping and filtering on subsets of intervals, these can be merged into a final VCF with the following.
#bcftools concat -n -o $workspace/final_merge.vcf.gz \
#	$workspace/interval1_final.vcf.gz \
#	$workspace/interval2_final.vcf.gz \
#	$workspace/interval3_final.vcf.gz 


#### Output individual (observed) heterozygosity.
#### This produces VCF files for each individual, then filters these based on missing data, star alleles, and minimum and maximum depth, then calls each confidently genotyped site as HomRef, HomAlt, or Het. 
#### This uses a default min depth of 15 and a max of 50. See Nielsen R et al 2011 Nat Rev Genet, and Li H 2014 Bioinformatics for details on setting these cutoffs. 
#### The second half of this code is built from user Kevin Blighe's (https://www.biostars.org/u/41557/) comment at https://www.biostars.org/p/298361/
#### This code can be easily extended for expected heterozygosity. See section below
for K in $samples;
do
	bcftools view -s ${K} $workspace/final.vcf.gz | \
		bcftools filter -e 'AVG(FMT/DP) < 15' | \
		bcftools filter -e 'ALT="*"' | \
		bcftools filter -e 'AVG(FMT/DP) > 50' | \
		bcftools norm -a --atom-overlaps . -Oz -o $workspace/${K}.vcf.gz 

	paste <(bcftools view $workspace/${K}.vcf.gz | \
		awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
		!/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
		\
	<(bcftools query -f '[\t%SAMPLE=%GT]\n' $workspace/${K}.vcf.gz | \
		awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}') \
		\
	<(bcftools query -f '[\t%SAMPLE=%GT]\n' $workspace/${K}.vcf.gz | \
		awk 'BEGIN {print "nHomAlt"} {print gsub(/1\|1|1\/1/, "")}') \
		\
	<(bcftools query -f '[\t%SAMPLE=%GT]\n' $workspace/${K}.vcf.gz | \
		awk 'BEGIN {print "nHomRef"} {print gsub(/0\|0|0\/0/, "")}') \
		\
	<(bcftools view $workspace/${K}.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
		!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
		\
	<(bcftools view $workspace/${K}.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
		!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\|1|1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
		\
	<(bcftools view $workspace/${K}.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
		!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|0|0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
		\
	| sed 's/,\t/\t/g' | sed 's/,$//g' > $workspace/${K}_het.txt
done



#### OPTIONAL: Output expected heterozygosity for a group of samples, "POP1".
#### This produces a VCF file for the group, then filters it based on missing data, spanning deletions, star alleles, and minimum and maximum depth, then for each confidently genotyped site lists the number of individuals that are HomRef, HomAlt, or Het. 
#### Sites with more than two alleles are atomised, retaining all variants for heterozygosity estimation.
#### Set the POP1_samples variable below before running this code

#POP1_samples= "Sample1,Sample2,Sample3,Sample4,Sample5"
#bcftools view -s $POP1_samples $workspace/final.vcf.gz | \
#	bcftools filter -e 'GT[*] = "mis"' | \
#	bcftools filter -e 'AVG(FMT/DP) < 15' | \
#	bcftools filter -e 'ALT="*"' | \
#	bcftools filter -e 'AVG(FMT/DP) > 50' | \
#	bcftools norm -a --atom-overlaps . -Oz -o $workspace/POP1.vcf.gz 
	

#paste <(bcftools view POP1.vcf.gz | \
#	awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
#	  !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
#	\
#  <(bcftools query -f '[\t%SAMPLE=%GT]\n' POP1.vcf.gz | \
#	awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}') \
#	\
#  <(bcftools query -f '[\t%SAMPLE=%GT]\n' POP1.vcf.gz | \
#	awk 'BEGIN {print "nHomAlt"} {print gsub(/1\|1|1\/1/, "")}') \
#	\
#  <(bcftools query -f '[\t%SAMPLE=%GT]\n' POP1.vcf.gz |\
#	awk 'BEGIN {print "nHomRef"} {print gsub(/0\|0|0\/0/, "")}') \
#	\
#  <(bcftools view POP1.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
#	!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
#	\
#  <(bcftools view POP1.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
#	!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\|1|1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
#	\
#  <(bcftools view POP1.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
#	!/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|0|0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
#	\
#  | sed 's/,\t/\t/g' | sed 's/,$//g' > $workspace/POP1_het.txt
