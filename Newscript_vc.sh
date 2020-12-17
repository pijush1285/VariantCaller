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


prefix=pat014_B_L3

#mkdir "$prefix"
#cd "$prefix"

Re1=/data/sata_data/workshop/wsu28/reference19/hg19.fa
r1=/data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs/pat014_B/pat014_B_TruSeq_AC4L8TACXX_L3_R1.fastq.gz
r2=/data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs/pat014_B/pat014_B_TruSeq_AC4L8TACXX_L3_R2.fastq.gz


#The alignment operation is performed using the code given below.
#ALignment
bwa mem -M -t 50 "$Re1" "$r1" "$r2" > "$prefix"_bwa.sam


#The mapped reads are now converted from sam file to bam file.
#Bam_File #Running
samtools view -S -b "$prefix"_bwa.sam > "$prefix"_bwa.bam -@ 50


#Now sort the Bam file using the code given below.
#Sort_Bam
samtools sort "$prefix"_bwa.bam > "$prefix"_Sort_bwa.bam


#Marked Duplicates
picard MarkDuplicates \
      I="$prefix"_Sort_bwa.bam \
      O="$prefix"_marked_duplicates.bam \
      M="$prefix"_marked_dup_metrics.txt


#Add replace Group
picard AddOrReplaceReadGroups \
       I="$prefix"_marked_duplicates.bam \
       O="$prefix"_RG_bwa.bam \
       RGID=4 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=20


#Creating the index file
picard BuildBamIndex I="$prefix"_RG_bwa.bam


#Variant call. Using the code given below you can call the variants.
# gatk --java-options "-Xmx10g" HaplotypeCaller \
#   -R "$Re1" \
#   -I "$prefix"_RG_bwa.bam \
#   -O "$prefix"_BWA_GATK.vcf.gz \
#   --native-pair-hmm-threads 50



# Or you can produce the GVCF
gatk --java-options "-Xmx4g" HaplotypeCaller \
  -R "$Re1" \
  -I "$prefix"_RG_bwa.bam \
  -O "$prefix".g.vcf.gz \
  --native-pair-hmm-threads 50 \
  -ERC GVCF




# https://gatk.broadinstitute.org/hc/en-us/articles/360036899732-GenotypeGVCFs
