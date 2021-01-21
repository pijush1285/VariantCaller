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


#Input path need to be edited
id=/data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs/pat154_B

#Output path need to be edited
op=/data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat154_B

count=0
for entry in "$id"/*
do
  #echo "$entry"
  array[$count]="$entry"
  count=$(($count+1))
  #echo $count
done

#Getting the number of samples in the Directory
counter=${#array[@]}
#counter=10

for ((number = 0; number < counter; number += 2));
do

    #For test dataset
    #word=$(awk -F/ '{print $8}' <<< "${array[$number]}")
    #echo $word
    #prefix=$(awk -F_ '{print $1}' <<< "${word}")_$(awk -F_ '{print $2}' <<< "${word}")_$(awk -F_ '{print $3}' <<< "${word}")
    #echo $prefix
    #For actual data set. I need to be edited. Please check before run
    word=$(awk -F/ '{print $9}' <<< "${array[$number]}")
    echo $word
    prefix=$(awk -F_ '{print $1}' <<< "${word}")_$(awk -F_ '{print $2}' <<< "${word}")_$(awk -F_ '{print $4}' <<< "${word}")_$(awk -F_ '{print $5}' <<< "${word}")
    echo $prefix

    r1="${array[$(($number))]}"
    echo The sample R2: $r1
    r2="${array[$(($number+1))]}"
    echo The sample R2: $r2

    #New output Paths
    mkdir "$op/$prefix"
    nop=$op/$prefix

    # Reference genome
    Re1=/data/sata_data/workshop/wsu28/reference19/hg19.fa


    #The alignment operation is performed using the code given below.
    #ALignment
    bwa mem -M -t 50 "$Re1" "$r1" "$r2" > "$nop/$prefix"_bwa.sam


    #The mapped reads are now converted from sam file to bam file.
    #Bam_File #Running
    samtools view -S -b "$nop/$prefix"_bwa.sam > "$nop/$prefix"_bwa.bam -@ 50


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
      --native-pair-hmm-threads 50 \
      -ERC GVCF



done




# Now for Combining all those VCF files the code is given below.

gatk CombineGVCFs \
-R /data/sata_data/workshop/wsu28/reference19/hg19.fa \
--variant /data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat014_B/pat014_B_L1/pat014_B_L8a.g.vcf.gz \
--variant /data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat014_B/pat014_B_L1a/pat014_B_L1a.g.vcf.gz \
--variant /data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat014_B/pat014_B_L2/pat014_B_L2.g.vcf.gz \
--variant /data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat014_B/pat014_B_L3/pat014_B_L3.g.vcf.gz \
--variant /data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat014_B/pat014_B_L4/pat014_B_L4.g.vcf.gz \
--variant /data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat014_B/pat014_B_L5/pat014_B_L5.g.vcf.gz \
--variant /data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat014_B/pat014_B_L6/pat014_B_L6.g.vcf.gz \
--variant /data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat014_B/pat014_B_L7/pat014_B_L7.g.vcf.gz \
--variant /data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat014_B/pat014_B_L8/pat014_B_L8.g.vcf.gz \
--variant /data/sata_data/workshop/wsu28/NidhanDaTask/human_wgs_cv/pat014_B/pat014_B_L8a/pat014_B_L1.g.vcf.gz \
-O Cohot_pat014_B.g.vcf.gz


# Now the gVCF need to be converted in to VCF files. The code is given below. 

gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /data/sata_data/workshop/wsu28/reference19/hg19.fa \
   -V Cohot_pat014_B.g.vcf.gz \
   -O Final_pat014.vcf.gz


