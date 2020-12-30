#!/bin/bash

################################################################################
#Now create the index of the reference file using bwa.
#Indexing
# bwa index "$Re1"

#Index file generate of the reference file.
# picard CreateSequenceDictionary \
#             R="$Re1" \
#             O=hg19.dict

# samtools faidx "$Re1"
################################################################################


for L in 1 2 3 4
do

  prefix=pat122_T_L"$L"



  #Paths
  ip=/data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs/pat122_T
  op=/data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat122_T
  mkdir "$op/$prefix"
  nop=$op/$prefix

  #Sample names example
  #r1=pat014_B_TruSeq_AC4L8TACXX_L"$L"_R1.fastq.gz
  #r2=pat014_B_TruSeq_AC4L8TACXX_L"$L"_R2.fastq.gz

  r1=pat122_T_TruSeqCD_BHW2LKDSXX_L"$L"_R1.fastq.gz
  r2=pat122_T_TruSeqCD_BHW2LKDSXX_L"$L"_R2.fastq.gz


  Re1=/data/sata_data/workshop/wsu28/reference19/hg19.fa


#The alignment operation is performed using the code given below.
#ALignment
bwa mem -M -t 25 "$Re1" "$ip/$r1" "$ip/$r2" > "$nop/$prefix"_bwa.sam


#The mapped reads are now converted from sam file to bam file.
#Bam_File #Running
samtools view -S -b "$nop/$prefix"_bwa.sam > "$nop/$prefix"_bwa.bam -@ 25


#Now sort the Bam file using the code given below.
#Sort_Bam
samtools sort "$nop/$prefix"_bwa.bam > "$nop/$prefix"_Sort_bwa.bam


#Marked Duplicates
picard MarkDuplicates \
      I="$nop/$prefix"_Sort_bwa.bam \
      O="$nop/$prefix"_marked_duplicates.bam \
      M="$nop/$prefix"_marked_dup_metrics.txt


#Add replace Group
picard AddOrReplaceReadGroups \
       I="$nop/$prefix"_marked_duplicates.bam \
       O="$nop/$prefix"_RG_bwa.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=20


#Creating the index file
picard BuildBamIndex I="$nop/$prefix"_RG_bwa.bam


#Variant call. Using the code given below you can call the variants.
# gatk --java-options "-Xmx10g" HaplotypeCaller \
#   -R "$Re1" \
#   -I "$prefix"_RG_bwa.bam \
#   -O "$prefix"_BWA_GATK.vcf.gz \
#   --native-pair-hmm-threads 50



# Or you can produce the GVCF
gatk --java-options "-Xmx4g" HaplotypeCaller \
  -R "$Re1" \
  -I "$nop/$prefix"_RG_bwa.bam \
  -O "$nop/$prefix".g.vcf.gz \
  --native-pair-hmm-threads 25 \
  -ERC GVCF


done

# https://gatk.broadinstitute.org/hc/en-us/articles/360036899732-GenotypeGVCFs
