## Overview
Completing functional annotation of selected *Brucella* isolates from Australian rodents. Isolates in question are:
- 0449, 0451, 0453, 0454, 0456, 0458, 0459, and 0460  

Previous work suggested that samples 0452, 0455, and 0457 are from an Acinetobacter sp., while 0450 is from a B. suis biovar 1 subtype. Functional analsis is restricted to those suspected to be related to *B. inopinata*.  

## Getting started: load the module and rename the files
Load the necessary module within Premise.  
```
module load linuxbrew/colsa
```

Fastq files were copied from the Cobb directory and reflect the generic naming structure applied by our sequencing center. Let's rename the fastq files to simplify basenames:
- wanted to rename them for less complicated basename... this cuts out all the index and L001/2 info:  
```
for i in *; do mv "$i" "`echo $i | sed "s/_.*L.*R/_R/"`"; done
```
- and this cuts out the unnecessary 001/002 info:  
```
for i in *; do mv "$i" "`echo $i | sed "s/_00*.//"`"; done
```
- and finally because we don't know what Brucella species to call them:  
```
rename Bsuis Brucella B*
```

We're going to trim the adapter sequences from these raw reads and apply some quality filtering next. Note that the parent directory `fastq` contains both the `raw` and `trim` subdirectories. To simplify the process of the forthcoming shell script I just made a symbolic link for each of the raw reads in the `trim` directory which will get deleted once the script is finished.  

## Adapter trimming with Skewer
To work with the SLURM job submission program we need to create two files: **slurm-skewer.sh** contains the commands needed to submit the job, while **skewer.sh** contains the skewer commands to actually perform the read trimming. Each are outlined separaetely below.  
First, we'll create a shell script which will execute the same skewere commands on all read pairs. Script is named `skewer.sh`.  
```
#!/bin/bash
cd /mnt/lustre/macmaneslab/devon/bruce/fastq/trim
ID_LIST=`ls -1 | sed 's/_.*//' | sort -u`
for READ in $ID_LIST
do
	skewer \
	--compress \
	--min 100 \
	--max 252 \
	--mode pe \
	--output ${READ} \
	--end-quality 10 \
	--mean-quality 5 \
	-x /mnt/lustre/macmaneslab/shared/adapters/TruSeq3-PE-2.fa \
	${READ}_R1.fastq.gz ${READ}_R2.fastq.gz
done
```

Next is the job submission script pointing to the directory of interest:
```
#!/bin/bash
#SBATCH -D /mnt/lustre/macmaneslab/devon/bruce/scripts
#SBATCH -p macmanes,shared
#SBATCH --job-name="oro-skewer"
#SBATCH --ntasks=1
#SBATCH --output=skewer.log

module purge
module load linuxbrew/colsa

srun skewer.sh
```

## Aligning reads
There are two steps when using the Burrows-Wheeler Aligner (BWA) program: indexing the reference genome, and then performing the alignment itself. Becuase you only index a single reference, I've put that command into a single slurm-executed shell script; however the alignment program is split between a job submission script and an alignment script. I've also created symbolic links within a new parent directory `bwa` to contain the trimmed .fq files.

Indexing the reference, *Brucella* sp. 83-13:  
```
#!/bin/bash
#SBATCH -D /mnt/lustre/macmaneslab/devon/bruce/bwa
#SBATCH -p macmanes,shared
#SBATCH --job-name="oro-bwaIndex"
#SBATCH --ntasks=1
#SBATCH --output=bwaIndex.log

module purge
module load linuxbrew/colsa

srun bwa index Brucella_sp_8313.fasta
```

Then creating a script for running the alignment program using the trimmed .fq files:  
```
#!/bin/bash
cd /mnt/lustre/macmaneslab/devon/bruce/bwa
ID_LIST=`ls *.gz | grep "trim" | sed 's/-.*//' | sort -u`
for READ in $ID_LIST
do
	bwa mem Brucella_sp_8313.fasta \
  -M \
  -t 24 \
  ${READ}-trimmed-pair1.fastq.gz ${READ}-trimmed-pair2.fastq.gz > \
  ${READ}.sam
done
```

The slurm submission script (not shown) followed just like the skewer submission script with the substitution of job and log file names.


