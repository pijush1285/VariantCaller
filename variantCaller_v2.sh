#!/bin/bash

echo -e "-----------------------------------------------------------------------------"
echo -e "                                                                             "
echo -e "This script will perform the variant calling operation step by step !"
echo -e "                                                                             "
echo -e "-----------------------------------------------------------------------------"



######################################################################################
#            This portion of code is used for alignment and fastQC
#####################################################################################
#Location of the tools
bwaMemLoc=/data/sata_data/workshop/wsu28/packages/bwa/bwa
fastqc=/data/sata_data/workshop/wsu28/packages/FastQC/fastqc
samtools=/data/sata_data/workshop/wsu28/packages/samtools1/samtools/samtools
bamtools=/data/sata_data/workshop/wsu28/anaconda3/envs/pdas/bin/bamtools
picard=/data/sata_data/workshop/wsu28/packages/picard/picard.jar
gatk=/data/sata_data/workshop/wsu28/packages/gatk/gatk


#Reference genome location
refFasta=/data/sata_data/workshop/wsu28/reference/hg38.fa


#For concatinating the fasta files used the code given below.
#cat U0a_CGATGT_L001_R1_00*.fastq > U0a_CGATGT_L001_R1_Marged.fastq
#cat U0a_CGATGT_L001_R2_00*.fastq > U0a_CGATGT_L001_R2_Marged.fastq

##############################################################################################
#                                 Location of the raw reads
##############################################################################################
longReadsL1R1S1=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R1_001.fastq
longReadsL1R1S2=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R1_002.fastq
longReadsL1R1S3=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R1_003.fastq
longReadsL1R1S4=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R1_004.fastq
longReadsL1R1S5=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R1_005.fastq

longReadsL1R2S1=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R2_001.fastq
longReadsL1R2S2=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R2_002.fastq
longReadsL1R2S3=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R2_003.fastq
longReadsL1R2S4=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R2_004.fastq
longReadsL1R2S5=/data/sata_data/workshop/wsu28/rawdata/Sample_U0a/U0a_CGATGT_L001_R2_005.fastq

echo "-------------------------------------------"
echo "All the sample files are taken as an Input"
echo "-------------------------------------------"
######################################################################
#                               FastQC
######################################################################
#The output file will be assembled under a directory
OutPutPath="/data/sata_data/workshop/wsu28/vcResults/"
samplelabel="Sample_U0aV2"

#Frist create Sample output folder
if [ ! -d $OutPutPath$samplelabel ]; then
  mkdir -p $OutPutPath$samplelabel;
fi

#Make a FastQc folder to save the FastQc results
FQoutDir="1_FastQC"
if [ ! -d $OutPutPath$samplelabel/$FQoutDir ]; then
  mkdir -p $OutPutPath$samplelabel/$FQoutDir;
fi

#cd /data/sata_data/workshop/wsu28/vcResults/FastQC
OutPutDir=$OutPutPath$samplelabel/$FQoutDir


logPaths=("$longReadsL1R1S1" "$longReadsL1R1S2" "$longReadsL1R1S3" "$longReadsL1R1S4" "$longReadsL1R1S5")
logPaths1=("$longReadsL1R2S1" "$longReadsL1R2S2" "$longReadsL1R2S3" "$longReadsL1R2S4" "$longReadsL1R2S5")


for i in ${!logPaths[@]};
do
    $fastqc ${logPaths[$i]} -o $OutPutDir
    $fastqc ${logPaths1[$i]} -o $OutPutDir
    echo "Processing complete: $i times"
done

echo "-------------------------------------------"
echo "FastQC operation has been compleated.."
echo "-------------------------------------------"

###################################################################################
#              Three type of indexing using different tools
###################################################################################

#1 BWA Indexing
#$bwaMemLoc index $refFasta


#2 Index by using samtools
#$samtools faidx $refFasta


#3 Index by using picard
#java -jar $picard CreateSequenceDictionary \
#REFERENCE=$refFasta \
#OUTPUT=hg38.dict


##################################################################################
#			        	Mapping
##################################################################################

#Make a MapFileContainer folder to save the mapped sam files results
MapFileContainer="2_MapFileContainer"
if [ ! -d $OutPutPath$samplelabel/$MapFileContainer ]; then
  mkdir -p $OutPutPath$samplelabel/$MapFileContainer;
fi

#Now change the directory.
OutPutDir1=$OutPutPath$samplelabel/$MapFileContainer



#Output File name
fileName=U0a_CGATGT_Default_marged

#Number of thread used
thread=20


#accuracy == "default"
cd $OutPutDir1
$bwaMemLoc mem $refFasta $longReadsL1R1S1 $longReadsL1R2S1 -t $thread > SplL1R1R2S1.out.sam
$bwaMemLoc mem $refFasta $longReadsL1R1S2 $longReadsL1R2S2 -t $thread > SplL1R1R2S2.out.sam
$bwaMemLoc mem $refFasta $longReadsL1R1S3 $longReadsL1R2S3 -t $thread > SplL1R1R2S3.out.sam
$bwaMemLoc mem $refFasta $longReadsL1R1S4 $longReadsL1R2S4 -t $thread > SplL1R1R2S4.out.sam
$bwaMemLoc mem $refFasta $longReadsL1R1S5 $longReadsL1R2S5 -t $thread > SplL1R1R2S5.out.sam


#accuracy == "strict"
#$bwaMemLoc mem -k 16 -W 60 $refFasta $longReads1 $longReads2 -t $thread > $fileName.out.sam

echo "-------------------------------------------"
echo "Mapping operation has been done......"
echo "-------------------------------------------"

############################################################################################
#                        AddOrReplaceReadGroups (Picard)
############################################################################################

logPaths2=("SplL1R1R2S1.out.sam" "SplL1R1R2S2.out.sam" "SplL1R1R2S3.out.sam" "SplL1R1R2S4.out.sam" "SplL1R1R2S5.out.sam")
logPaths3=("SplL1R1R2S1.bam" "SplL1R1R2S2.bam" "SplL1R1R2S3.bam" "SplL1R1R2S4.bam" "SplL1R1R2S5.bam")

for i in ${!logPaths2[@]};
  do
    #Input file
    input=$OutPutDir1/${logPaths2[$i]}
    echo "----------------------------------"
    echo $input
    echo "Input file name and location"
    #Output file
    output=$OutPutDir1/${logPaths3[$i]}
    echo "----------------------------------"
    echo $output
    echo "Output file name and location"

    java -jar $picard AddOrReplaceReadGroups \
          I=$input \
          O=$output \
          RGID=4 \
          RGLB=lib1 \
          RGPL=ILLUMINA \
          RGPU=unit1 \
          RGSM=20

          echo "Processing complete: $i times"

done

echo "----------------------------------------------------"
echo "Add remove group operation has been done ..         "
echo "----------------------------------------------------"

#############################################################################################
#                 Base Quality Score Recalibration (Create table)
#############################################################################################

#Make a BQSRcreateTable folder to save the mapped sam files results
BQSRcreateTable="3_BQSRcreateTable"
if [ ! -d $OutPutPath$samplelabel/$BQSRcreateTable ]; then
  mkdir -p $OutPutPath$samplelabel/$BQSRcreateTable;
fi

#New Directory location
OutPutDir2=$OutPutPath$samplelabel/$BQSRcreateTable


#Before using the vcf file you would required to make a tab index vcf file
#tabix -p vcf 1000G_phase1.snps.high_confidence.hg38.vcf.gz

goldRndelsVCF=/data/sata_data/workshop/wsu28/vcf1/1000G_phase1.snps.high_confidence.hg38.vcf.gz
setOfSitesToMask=/data/sata_data/workshop/wsu28/vcf1/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

#Input file names for this function
logPath4=("SplL1R1R2S1.bam" "SplL1R1R2S2.bam" "SplL1R1R2S3.bam" "SplL1R1R2S4.bam" "SplL1R1R2S5.bam")
logPath5=("SplL1R1R2S1" "SplL1R1R2S2" "SplL1R1R2S3" "SplL1R1R2S4" "SplL1R1R2S5")



for i in ${!logPath4[@]};
  do

    #Location of input files (This is the input file)
    dedupReads=$OutPutDir1/${logPath4[$i]}
    echo "-----------------------------------------------------------------"
    echo ${logPath4[$i]}
    echo "BAM file reading compleated ......"
    echo "-----------------------------------------------------------------"


    #Output File name
    tgList=$OutPutDir2/${logPath5[$i]}-recal_data.table.txt
    tgList1=$OutPutDir2/${logPath5[$i]}-recalibrate.out
    tgList2=$OutPutDir2/${logPath5[$i]}-recalibrate.err


    #Now using the base recalibration
    $gatk BaseRecalibrator \
      -I $dedupReads \
      -R $refFasta \
      --known-sites $goldRndelsVCF \
      --known-sites $setOfSitesToMask \
      -O $tgList \
      >$tgList1 2>$tgList2

      echo "-----------------------------------------------------------------"
      echo "Base Quality Score Recalibration operation has been done ......"
      echo "-----------------------------------------------------------------"

done


#############################################################################################
#                    Base Quality Score Recalibration (apply)
#############################################################################################

#Make a BQSRcreateTable folder to save the mapped sam files results
BQSRcreateApply="4_BQSRcreateApply"
if [ ! -d $OutPutPath$samplelabel/$BQSRcreateApply ]; then
  mkdir -p $OutPutPath$samplelabel/$BQSRcreateApply;
fi

#New Directory location
OutPutDir3=$OutPutPath$samplelabel/$BQSRcreateApply


#Input file names for this function
logPath4=("SplL1R1R2S1.bam" "SplL1R1R2S2.bam" "SplL1R1R2S3.bam" "SplL1R1R2S4.bam" "SplL1R1R2S5.bam")
logPath5=("SplL1R1R2S1" "SplL1R1R2S2" "SplL1R1R2S3" "SplL1R1R2S4" "SplL1R1R2S5")


for i in ${!logPath4[@]};
  do

    #Location of input files (This is the input file)
    dedupReads=$OutPutDir1/${logPath4[$i]}
    tgList=$OutPutDir2/${logPath5[$i]}-recal_data.table.txt

    echo "-----------------------------------------------------------------"
    echo ${logPath4[$i]}
    echo "BAM file reading compleated ......"
    echo "-----------------------------------------------------------------"


    #Output File name

    tgLista=$OutPutDir3/${logPath5[$i]}-GATK.bam
    tgList1a=$OutPutDir3/${logPath5[$i]}-recalbam.out
    tgList2a=$OutPutDir3/${logPath5[$i]}-recalbam.err


    $gatk PrintReads \
        -R $refFasta \
        -I $dedupReads \
        -O $tgLista \
        >$tgList1a 2>$tgList2a


        echo "---------------------------------------------------------------------------"
        echo "Base Quality Score Recalibration (apply) operation has been done .."
        echo "---------------------------------------------------------------------------"

done


#############################################################################################
#                                Now gather all the BAM files.
#############################################################################################

cd $OutPutDir3
cat SplL1R1R2S*-GATK.bam > GatherAllBAMs.bam

echo "---------------------------------------------------------------------------"
echo "The BAM files have been gather..."
echo "---------------------------------------------------------------------------"

#############################################################################################
#                       Now break the GatherAllBAMs.bam file chromosom wise.
#############################################################################################

#Now make a folder which contains the broken bam file chromosom wise.
ChromosomWiseBreak="5_ChromosomWiseBreak"
if [ ! -d $OutPutPath$samplelabel/$ChromosomWiseBreak ]; then
  mkdir -p $OutPutPath$samplelabel/$ChromosomWiseBreak;
fi

#New Directory location
OutPutDir4=$OutPutPath$samplelabel/$ChromosomWiseBreak

#Copy the GatherAllBAMs.bam file in to ChromosomWiseBreak directory
cp $OutPutDir3/GatherAllBAMs.bam $OutPutDir4

#Now change the directory
cd $OutPutDir4

#Now break the bam file chromosom wise
$bamtools split -in GatherAllBAMs.bam -reference

echo "---------------------------------------------------------------------------"
echo "Now break the GatherAllBAMs.bam file chromosom wise operation is done .."
echo "---------------------------------------------------------------------------"

#Now remove the GatherAllBAMs file from this folder.
rm GatherAllBAMs.bam


############################################################################################
#              Create bam file index for the GATK bam file
############################################################################################


#Now make a folder which contains the broken bam file chromosom wise.
GatherBam="5_GatherBam"
if [ ! -d $OutPutPath$samplelabel/$GatherBam ]; then
  mkdir -p $OutPutPath$samplelabel/$GatherBam;
fi


#New Directory location
OutPutDir5=$OutPutPath$samplelabel/$GatherBam

#Copy the GatherAllBAMs.bam file in to ChromosomWiseBreak directory
cp $OutPutDir3/GatherAllBAMs.bam $OutPutDir5

#Now change the directory
cd $OutPutDir5



#Sorted the align sam file to converted into sorted reads
java -jar $picard SortSam \
     I=GatherAllBAMs.bam \
     O=GatherAllBAMs_sorted.bam \
     SORT_ORDER=coordinate


#Now creating the index file.
java -jar $picard BuildBamIndex \
      -I GatherAllBAMs_sorted.bam \
      -O GatherAllBAMs.bai \
      >GatherAllBAMs_bamindex.out 2>GatherAllBAMs_bamindex.err


echo "BAM file indexing operation is ended......"
echo "--------------------------------------------------"


############################################################################################
#              Now Haplotype call operation using GATK
############################################################################################

#Now Haplotype call operation will be performed.
$gatk --java-options "-Xmx4g" HaplotypeCaller \
     -R $refFasta \
     -I GatherAllBAMs.bam \
     --dbsnp $goldRndelsVCF \
     -O GatherAllBAMs.vcf.gz \
     -bamout GatherAllBAMs_bamout.bam



echo "HaplotypeCaller calling operation is ended....."
echo "--------------------------------------------------"

#############################################################################################
