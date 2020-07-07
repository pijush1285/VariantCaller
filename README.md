# Variant-calling-for-the-beginners

**Variant Calling:**
Variant calling is the mechanism through which variants from the sequence data are found. This process can be exicuted in to three steps. Those steps are given below.
+ Execute entire genome or exome sequences to create FASTQ files.
+ Align the sequences by generating BAM or CRAM files to a reference genome.
+ Identify and write to a VCF file where the aligned reads differ from the reference genome.

This pipeline will show you how to make a variant call step by step using different publickly available tools.
Initially a .sh file will be given which contains a newly build variant calling pipeline.
Gradually the whole section will be illustrated.


<font size="5">**1. TOOLS**</font>

There are lots of available tools for varient calling. Some combination of tools produces best variant call result.
In this section I enlist the tools which are used for variant calling operation.

+ **GATK** (*Download Url: https://github.com/broadinstitute/gatk*)
+ **picard** (*Download Url: https://github.com/broadinstitute/picard*)
+ **samtools** (*Download Url: https://github.com/samtools/samtools*)
+ **bwa** (*Download Url: https://github.com/lh3/bwa*)
+ **bamtools** (*Download Url: https://github.com/pezmaster31/bamtools*)
+ **FastQC** (*Download Url: https://github.com/s-andrews/FastQC*)
+ **MultiQC** (*Download Url: https://github.com/ewels/MultiQC*)


<font size="5">**2. DATA SETS**</font>

+ *Raw Data sets*

A slanderous practice data set is provided here to run the entire variance calling the pipeline. The user can easily download the entire data set on their local server and run the pipeline easily.

 **Sample:** *Sample_U0a*
 **Project:** *Project_RM8398*
 **Sequencer:** *HiSeq_300x*
 **Download Url:** 

<font size="2">ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/Project_RM8398/Sample_U0a/</font>


+ *Reference data set*

A reference genome (also known as a reference assembly) is a database of the digital nucleic acid sequence, assembled by scientists as a representative example of a set of genes for a species. In this pipeline I have used GRCh38 for mapping the fastq files.

Download Url: https://www.ncbi.nlm.nih.gov/grc/human


+ *SNP data set*

Another essential thing that is needed for variant calling is the SNP database. A download link for downloading the SNP database is provided here below. 

Download Url: ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/


<font size="4">**3. Basic workflow of this pipeline**</font>

This pipeline works very simply, executing various functions from different packages, step by step. The novices just need to understand the functionality of those functions. After a few successful executions of this pipeline with different data sets, it will be more flexible for the user to control input arguments that are needed for these functions. Initially the user only require to change tha path of those data sets and run the pipeline using the command given below.

+ 1. Take a look at the FASTQ file

It is required to put a glance in the Fastq files. After that continue the processing. 

~~~
$./variantCaller_v2.sh
~~~


<font size="5">**4. Other Pipelines for variant call**</font>

Here, two more pipelines are included, which are very simple pipelines for calling variants.
The name of those pipelines are 1) BWA_samtools_bcftools.sh and 2) BWA_Gatk_v2.sh.


In most cases, these packages can be downloaded and built using this pattern:

~~~
git clone https://github.com/lh3/bwa.git 
cd bwa; make
~~~

<img src="Images/pipeline1.png">
