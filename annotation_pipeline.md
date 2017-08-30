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

