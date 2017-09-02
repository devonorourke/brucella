## getting started (in a perfect world)
`nasp` is available on Premise within the `anaconda/colsa` virtual environment, and you can load that collection of programs and scripts once you've logged in:
```
module purge
module load anaconda/colsa
source activate nasp-1.0.2
```
The program begins by just entering the term `nasp`, and will create a series of further command line prompts. So to begin:
```
nasp
```
Except this didn't work. The process kept crashing at the .bam --> .vcf step. Because the world ain't perfect.  

## getting started again (in an imperfect world)
Because I was having a lot of trouble with a few of the programs functioning, I wanted to edit one of the big python scripts (`dispatcher.py`). This can't be done within the default `nasp` virtual environment without admin privileges, but there is a nice workaround: just making a copy of the entire virtual environment and then editing whatever you want therein. Because this anaconda install is using Python3, we can do this by entering just a single line command and have our own personal **NASP** virtual environment and edit whatever scripts we want. See [this site](https://conda.io/docs/user-guide/tasks/manage-environments.html#cloning-an-environment) for further instructions. I executed this command in the directory within the generic `nasp` virtual environment (`$HOME/.conda/envs`).  
```
conda create --name mynasp --clone nasp-1.0.2
```

Once that virtual environment is installed you're going to want to edit the `dispatcher.py` located at `$HOME/.conda/envs/mynasp/lib/python3.6/site-packages/nasp`. Specifically, there are a pair of lines (lines 443 and 444) which should be substituted as follows:
```
## original code to change
command_parts = ["%s mpileup -uD -d 10000000 -f %s %s" % (sampath, reference, bam_file),
                     "%s view -ceg %s - > %s" % (path, args, final_file)]

## edited code to insert
command_parts = ["%s mpileup -u -t DP -d 8000 -f %s %s" % (sampath, reference, bam_file),
                     "bcftools call -c > %s" % (final_file)]
```
It turns out we were using an unnecessary depth of reads in the pileup file (`-d 1000000` is ridiculous). Moreover, we were using an older series of flags within the `samtools view` argument to create the .vcf file call, but those don't exist anymore with the update samtools package that is installed within this `nasp` virtual environment package. Once those edits are made, save that file (overwrite the existing one), and you'll be set to start running the program.

## On to actually running NASP!

Thereafter you need to simply execute a series of filepaths (always use the full path) as well as a few other questions. In general we're going to use the default questions through these prompts. All responses for the various `nasp` runs are indicated below. 

### Broad analysis:
These commands were used for the "broad" output - using a handful of Brucella representing many disparate clades. Followed default setting when not specified. 
```
#output:
/mnt/lustre/macmaneslab/devon/bruce/nasp/results_all/broad_results

#reference path:
/mnt/lustre/macmaneslab/devon/bruce/nasp/rodent_ref/Brucella_sp_8313.fasta

#job manager?
SLURM
  #partition?
  macmanes
 
#external references:
/mnt/lustre/macmaneslab/devon/bruce/nasp/external_refs_all/broad_refs

#read files?
N

#pre-called VCFs?
Y
  #where?
  /mnt/lustre/macmaneslab/devon/bruce/nasp/results_all/narrow_results/samtools
```

### Narrow analysis
These settings were used for the narrow experiment: narrow in the sense that we only provided a single external genome (Brucella_sp_NF2653) instead of the many various Brucella samples. Same reference used for alignment (B.sp_8313). Followed default setting when not specified. 

```
#output:
/mnt/lustre/macmaneslab/devon/bruce/nasp/results_all/narrow_results

#reference path:
/mnt/lustre/macmaneslab/devon/bruce/nasp/rodent_ref/Brucella_sp_8313.fasta

#job manager?
SLURM
  #partition?
  macmanes

#external references:
/mnt/lustre/macmaneslab/devon/bruce/nasp/external_refs_all/narrow_refs

#read fastq files:
Y

#where?
/mnt/lustre/macmaneslab/devon/bruce/nasp/australia_fq

#run Bowtie2?
N

#run GATK?
N

#run VarScan?
N
```


nasp-1.0.2/lib/python3.6/site-packages/nasp/nasptool_linux_64
