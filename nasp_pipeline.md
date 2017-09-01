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

Once that virtual environment is installed you're going to want to edit the `dispatcher.py` located at `$HOME/.conda/envs/mynasp/lib/python3.6/site-packages/nasp` with the following things in mind:

## On to actually running NASP!

Thereafter you need to simply execute a series of filepaths (always use the full path) as well as a few other questions. In general we're going to use the default questions through these prompts. All responses for the various `nasp` runs are indicated below. 

### Broad analysis:
These commands were used for the "broad" output - using a handful of Brucella representing many disparate clades.
```
#output:
/mnt/lustre/macmaneslab/devon/bruce/nasp/results_all/broad_results

#reference path:
/mnt/lustre/macmaneslab/devon/bruce/nasp/rodent_ref/Brucella_sp_8313.fasta

#check for dups?
Y

#job manager?
SLURM
  #partition?
  macmanes
  #additional requirements?
  {blank}

#external references:
/mnt/lustre/macmaneslab/devon/bruce/nasp/external_refs_all/broad_refs
  #advanced NUCmer?
  no

#read fasta files:
/mnt/lustre/macmaneslab/devon/bruce/nasp/australia_samps

#pre-aligned SAM/BAM files to include?
N

#pre-called VCFs?
N

#advanced matrix generator settings?
N
```

### Narrow analysis
These settings were used for the narrow experiment: narrow in the sense that we only provided a single external genome (Brucella_sp_NF2653) instead of the many various Brucella samples. Same reference used for alignment (B.sp_8313)

```
#output:
/mnt/lustre/macmaneslab/devon/bruce/nasp/results_all/narrow_results

#reference path:
/mnt/lustre/macmaneslab/devon/bruce/nasp/rodent_ref/Brucella_sp_8313.fasta

#check for dups and skip SNPs?
Y

#job manager?
SLURM
  #partition?
  macmanes
  #additional requirements?
  {default}

#external references:
/mnt/lustre/macmaneslab/devon/bruce/nasp/external_refs_all/narrow_refs
  #advanced NUCmer?
  no

#read fastq files:
N

#pre-aligned SAM/BAM files to include?
Y

#where?
/mnt/lustre/macmaneslab/devon/bruce/nasp/australia_bams

#pre-called VCFs?
N
  #min coverage threshold [10]?
  {default}
  #min acceptable proportion[0.9]?
  {default}

#advanced matrix generator settings?
N
```


nasp-1.0.2/lib/python3.6/site-packages/nasp/nasptool_linux_64
