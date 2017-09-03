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
There are two steps when using the [Burrows-Wheeler Aligner (BWA) program](http://bio-bwa.sourceforge.net/bwa.shtml): indexing the reference genome, and then performing the alignment itself. Becuase you only index a single reference, I've put that command into a single slurm-executed shell script; however the alignment program is split between a job submission script and an alignment script. I've also created symbolic links within a new parent directory `bwa` to contain the trimmed .fq files. **Note that these sym. links were manually removed following the completion of the script and won't appear in the directory**

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

## Converting .sam files to indexed .bam files
Side note - see [this page](https://davetang.org/wiki/tiki-index.php?page=SAMTools) for a million helpful one-liners when working with .sam and .bam files (and more).  
This process is a helpful step to save time and memory for generating a consensus fasta sequence. Besides the generic slurm script (not shown), there is a two-part shells script which first compresses the .sam to .bam file, then indexes that .bam file.  
```
#!/bin/bash
cd /mnt/lustre/macmaneslab/devon/bruce/bwa
ID_LIST=`ls *.sam | sed 's/.sam//' | sort -u`
for READ in $ID_LIST
do
  samtools view -bS ${READ}.sam | samtools sort - -o ${READ}.bam
done

for READ in $ID_LIST
do
  samtools index ${READ}.bam ${READ}.bai
done
```

**Note that the .sam files were removed following the creation of the indexed .bam and .bai files**. You can always convert .bam back to .sam.  

## Creating fasta files from .sam files (generating a consensus sequence):
See [this page](http://www.metagenomics.wiki/tools/samtools/consensus-sequence) for further info on where these commands were taken from. I had to install bcftools directly (not on Premise) which amounted to following [these directions](http://www.htslib.org/download/). I used the file existing directory as the location to install the `bin` folder containing the programs, then made symbolic links to my `$HOME/bin` directory for all those programs so that they could be called from my existing $PATH variable (ie. I didn't export any new $PATH).  
We're going to derive a consensus .fasta sequence by first converting our .bam to .fastq, then use the quality information in that .fq file to filter our .fa file a bit. Rather than splitting up the script into multiple sections we're going to just pipe everything in one long command. As with the last few commands, I'm not showing the slurm submission shell script.
```
#!/bin/bash
cd /mnt/lustre/macmaneslab/devon/bruce/bwa
ID_LIST=`ls *.bam | sed 's/.bam//' | sort -u`
for READ in $ID_LIST
do
samtools mpileup -uf Brucella_sp_8313.fasta ${READ}.bam | bcftools call -c | vcfutils.pl vcf2fq | \
seqtk seq -q20 -n N -A > /mnt/lustre/macmaneslab/devon/bruce/samp_fastas/${READ}.fa
done
```

## detour before Prokka --> making a .vcf file and annotating that instead of making Prokka do the work
Needed to create a .vcf file which basically takes a portion of the script above (calling a consensus fasta) and just cuts it at the `bcftools call ...` step, rather than taking that .vcf file and converting into a .fq/.fa file. First part of the detour involves generating that .vcf file; second part of the detour invovlves annotating that .vcf file directly with a package called [snpEff](http://snpeff.sourceforge.net/index.html).  

First, the .vcf file. Note that within the `bcftools call` argument we're specifying that we're dealing with haploids here.
```
#!/bin/bash
#SBATCH -D /mnt/lustre/macmaneslab/devon/bruce/scripts
#SBATCH -p macmanes,shared
#SBATCH --job-name="oro-vcfs"
#SBATCH --ntasks=1
#SBATCH --output=australia_vcfs.log

module purge
module load linuxbrew/colsa

cd /mnt/lustre/macmaneslab/devon/bruce/bwa

ID_LIST=`ls *.bam`

srun samtools mpileup -uf Brucella_sp_8313.fasta $ID_LIST | \
bcftools call --consensus-caller --variants-only --ploidy 1 > \
/mnt/lustre/macmaneslab/devon/bruce/vcfs/australia.vcf
```

Next, we're going to install snpEff and use it to annotate the .vcf file. See [install instructions](http://snpeff.sourceforge.net/download.html) for more details. Briefly:
```
## Install package and create soft link to .jar files in $HOME/bin
cd $HOME/pkgs
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
ln -s /mnt/lustre/macmaneslab/devon/pkgs/snpEff/*.jar /mnt/lustre/macmaneslab/devon/bin/
	## note you'll still need to specify the path to the .jar files as $HOME/bin (it just shortens it up a bit)
```

Create the database you want. First, download a list of databases available. Then perform a grep search to find which one you want. Then download it.
```
## download the list:
cd $HOME/bruce
mkdir snpEff
cd snpEff
java -jar ~/bin/snpEff.jar databases > databaseList.txt

## find the 83-13 reference Brucella:
grep -c "Brucella" databaseList.txt	#yields 332 entries
grep -c "Brucella_sp" databaseList.txt	#yields 16 entries
grep "Brucella_sp" databaseList.txt 	#our species is listed as "Brucella_sp_83_13"

## download your database:
java -jar ~/bin/snpEff.jar download Brucella_sp_83_13
```

Now we can annotate the .vcf file. See their [manual page](http://snpeff.sourceforge.net/SnpEff_manual.html) for loads of helpful flags to pass. I've created a short shell script to exeucte this with SLURM, but for such a tiny file you could really just do this from the head node directly on our compute cluster. Unlike other shell scripts, be sure to execute the script within `/mnt/lustre/macmaneslab/devon/bruce/snpEff/`, and specify the full path to the shell script (`/mnt/lustre/macmaneslab/devon/bruce/scripts/snpEff.sh`). In other words, create this shell file:

```
#!/bin/bash
#SBATCH -D /mnt/lustre/macmaneslab/devon/bruce/snpEff/
#SBATCH -p macmanes,shared
#SBATCH --job-name="oro-snpEff"
#SBATCH --ntasks=1
#SBATCH --output=snpEff.log

module purge
module load linuxbrew/colsa

cd /mnt/lustre/macmaneslab/devon/bruce/snpEff

srun java -jar ~/bin/snpEff.jar \
Brucella_sp_83_13 \
/mnt/lustre/macmaneslab/devon/bruce/vcfs/australia.vcf > \
/mnt/lustre/macmaneslab/devon/bruce/snpEff/australia_ann.vcf
```

Make the script executable (`chmod +x`) and then execute the command as follows:
```
cd /mnt/lustre/macmaneslab/devon/bruce/snpEff
sbatch /mnt/lustre/macmaneslab/devon/bruce/scripts/snpEff.sh
```
The output of this file provides the raw variant calls - note there are no depth, quality, or variant type filtered out among these calls - those are done next.  

### Filtering the .vcf file
We can apply a suite of filters to this annotated variant call file (`{name}.ann.vcf`). This is done using a related program to `snpEff` called **SnpSift**. This should be executed using a shell script for larger files but takes about 5 seconds for .vcf files that have less than 200 SNPs total... so what I've put below is the non-shell script portion of that code.  

#### Basic quality and depth filters
First, let's filter any sites which contain an average Phred-scaled quality score of less than 20 AND contain a depth of less than 100 (we usually have really high coverage for all sites, and this `100` is an aggregate of the 8 samples, so this is really more of an average of about 12 reads per site assuming that all 8 samples are present in the SNP).  
```
java -jar ~/bin/SnpSift.jar filter "( DP > 100 ) & ( QUAL >= 20 )" \
/mnt/lustre/macmaneslab/devon/bruce/snpEff/australia_ann.vcf \
> /mnt/lustre/macmaneslab/devon/bruce/snpEff/australia-qualfilt.ann.vcf
```
This results in a file containing 255 raw SNPs filtered down to 238 SNPs.  

#### Non synonymous sites only
There are a many different variables you could filter from the annotated .vcf file. For example, you may decide to filter by the type of variant described. To get a list of all the possible strings (names) that this variable could encompass you'll need to apply the little one-liner below. One thing to manipulate would be the second `cut` command, moving the `f` field delimiter around for the various annotations. In this instance, we're targeting the second chunk of the last field of the .vcf file (the annotation field), which describes whether the mutation is a missense mutation, synonymous variant, or something else:
```
grep -v "^#" australia-qualfilt.ann.vcf | awk '{print $8}' | cut -f 17 -d ';' | cut -f 2 -d '|' | sort -u
```
So we find there are a series of mutation types:  
- conservative_inframe_deletion
- disruptive_inframe_insertion
- frameshift_variant
- frameshift_variant&stop_gained
- missense_variant
- stop_gained
- stop_lost&splice_region_variant&conservative_inframe_deletion
- synonymous_variant
- upstream_gene_variant  

Two possible things to filter here:  
First, we might want to just keep the `missense_variant` SNPs *and nothing else*:
```
java -jar ~/bin/SnpSift.jar filter "ANN[*].EFFECT has 'missense_variant'" australia-qualfilt.ann.vcf > \
australia-missenseOnly.ann.vcf
```

Second, we might want to just filter **out** one specific variant type, say `synonymous_variant`, but *nothing else*:



## Annotation using Prokka
One word of caution... The output .fa files (for each Australian sample) have pipe's separating out the header elements - that is, they exactly match the ** *Brucella* sp. 83-13** input sample we provided. This shouldn't be a problem with running Prokka, but if you need to manipulate headers you can do that with a one-liner `sed` command as follows. Just make sure you apply it to **both the samples and the reference** fastas:

```
## changing headers in fasta (no specific file name given):
sed -i '/^>sid/s/ /-/g' {whatevername}.fa
```

We'll use Prokka as our gene annotation program for good reasons (because it's what's available!). See [here](https://github.com/tseemann/prokka#invoking-prokka) for some example commands from the program's Github page.  

```
#!/bin/bash
cd /mnt/lustre/macmaneslab/devon/bruce/samp_fasta
ID_LIST=`ls *.fasta | sed 's/.fasta//' | sort -u`
for READ in $ID_LIST
do
	prokka \
	--prefix ${READ} \
	--genus Brucella \
	--usegenus \
	--gram neg \
	--cpus 0 \
	{READ}.fasta
done
```

