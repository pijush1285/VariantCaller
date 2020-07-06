#!/bin/bash

echo -e "-----------------------------------------------------------------------------"
echo -e "                                                                             "
echo -e "This script will perform the variant calling operation step by step !"
echo -e "This program is developed for the analysis of the lab made dataset....       "
echo -e "                                                                             "
echo -e "-----------------------------------------------------------------------------"


# The code given below will ask the names of the reference file, read 1 and read 2 files.
read -p " enter your Reference:" Re1
read -p " enter the Read 1:" r1
read -p " enter the Read 2:" r2

#Or the user can directly assigne the files names into the variables
Re1=finalRefGen.fa
r1=NRead3.fq
r2=NRead4.fq

echo "-------------------------------------------"
echo "All the sample files are taken as an Input"
echo "-------------------------------------------"



#Location of the tools. Those location of the tools provided here is situated in our server.
#The user must have to insert the location of the tools where the packages are installed.
#You may use conda, in that case you do not have to specify where the tools are located.
bwa=/data/sata_data/workshop/wsu28/packages/bwa/bwa
samtools=/data/sata_data/workshop/wsu28/packages/samtools1/samtools/samtools
bcftools=/data/sata_data/workshop/wsu28/packages/samtools1/bcftools/bcftools
picard=/data/sata_data/workshop/wsu28/packages/picard/build/libs/picard.jar
gatk=/data/sata_data/workshop/wsu28/packages/gatk/gatk




#Now create the index of the reference file using bwa. 
#Indexing
$bwa index "$Re1"

echo "-------------------------------------------"
echo "Index of the reference file is created"
echo "-------------------------------------------"


#The alignment operation is performed using the code given below.
#ALignment 
$bwa mem -M -t 20 "$Re1" "$r1" "$r2" > "$r1"_"$r2"_bwa.sam

echo "-------------------------------------------"
echo "The reads are mapped with the reference genome."
echo "-------------------------------------------"


#The mapped reads are now converted from sam file to bam file.
#Bam_File #Running
$samtools view -S -b "$r1"_"$r2"_bwa.sam > "$r1"_"$r2"_bwa.bam -@ 20

echo "-------------------------------------------"
echo "Sam file to Bam file conversion is done."
echo "-------------------------------------------"


#Now sort the Bam file using the code given below.
#Sort_Bam 
$samtools sort "$r1"_"$r2"_bwa.bam > "$r1"_"$r2"_Sort_bwa.bam

echo "-------------------------------------------"
echo "Sorting operation is performed using samtools."
echo "-------------------------------------------"


#Marked Duplicates  
java -jar $picard MarkDuplicates \
      I="$r1"_"$r2"_Sort_bwa.bam \
      O="$r1"_"$r2"_marked_duplicates.bam \
      M="$r1"_"$r2"_marked_dup_metrics.txt
	  
echo "-------------------------------------------"
echo "Marked duplicate operation is performed"
echo "-------------------------------------------"


#Add replace Group   
java -jar $picard AddOrReplaceReadGroups \
       I="$r1"_"$r2"_marked_duplicates.bam \
       O="$r1"_"$r2"_RG_bwa.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=20
	   
echo "-------------------------------------------"
echo "Add replace group operation is performed"
echo "-------------------------------------------"



#Creating the index file   
jAdd replace Groupava -jar $picard BuildBamIndex \
      I="$r1"_"$r2"_RG_bwa.bam
	  
echo "-------------------------------------------"
echo "Index file is created"
echo "-------------------------------------------"




#Index file generate of the reference file.
java -jar $picard CreateSequenceDictionary \ 
      R="$Re1" \ 
      O=finalRefGen.dict
	  
echo "-------------------------------------------"
echo "Index file of the reference file is done"
echo "-------------------------------------------"



#Variant call. Using the code given below you can call the variants.
$gatk --java-options "-Xmx10g" HaplotypeCaller -R "$Re1" -I "$r1"_"$r2"_RG_bwa.bam -O "$r1"_BWA_GATK_"$r2".vcf.gz

echo "-------------------------------------------"
echo "Variant call operation is performed"
echo "-------------------------------------------"


#In case if you have more than one vcf files, you need to merge those file together. 
#The code given below only join together all the vcf files. After joining those file you will
#require to remove the duplicates.
#Merge VCF files. This code is working here.
java -jar $picard MergeVcfs I=NRead1.fq_BWA_GATK_NRead2.fq.vcf.gz I=NRead3.fq_BWA_GATK_NRead4.fq.vcf.gz O=Final_merged_variants.vcf.gz

#The code removes the duplicates.
bcftools norm -D Final_merged_variants.vcf.gz -o Final_merged_variants1.vcf.gz

echo "-------------------------------------------"
echo "Merged the vcf files and removing the duplicates"
echo "-------------------------------------------"



#the code given below can be used for Indel detection.
#variant_Sepration_Indel_SNV
#vcftools --vcf output/"$r1"/BWA_GATK_"$r1".vcf --remove-indels --recode --recode-INFO-all --out output/"$r1"/BWA_GATK_SNP_"$r1".vcf
#vcftools --vcf output/"$r1"/BWA_GATK_"$r1".vcf --keep-only-indels  --recode --recode-INFO-all --out output/"$r1"/BWA_GATK_Indels_"$r1".vcf



#Few code given below for compering vcf files.
#Total no of SNP
bcftools view -v snps Final_merged_variants1.vcf | grep -v "^#" | wc -l
bcftools view -v snps NRead1.fq_BWA_GATK_NRead2.fq.vcf.gz | grep -v "^#" | wc -l
bcftools view -v snps NRead3.fq_BWA_GATK_NRead4.fq.vcf.gz | grep -v "^#" | wc -l



#This code is also required.
bcftools view -v snps Final_merged_variants1.vcf | perl -lane 'if (/^#/) { print } elsif 
(length($F[3]) == 1) { print }' | bgzip > Final_merged_variants1.vcf.gz


#Index
tabix -p vcf Final_merged_variants1.vcf.gz

#Calculate Jaccard index
bedtools jaccard -a snpOriginal.vcf.gz -b Final_merged_variants1.vcf.gz



