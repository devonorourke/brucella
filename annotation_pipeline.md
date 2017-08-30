## Overview
Completing functional annotation of selected *Brucella* isolates from Australian rodents. Isolates in question are:
- 0449, 0451, 0453, 0454, 0456, 0458, 0459, and 0460  

Previous work suggested that samples 0452, 0455, and 0457 are from an Acinetobacter sp., while 0450 is from a B. suis biovar 1 subtype. Functional analsis is restricted to those suspected to be related to *B. inopinata*.  

## Getting started
Load the necessary module within Premise.  
```
module load linuxbrew/colsa
```

Now run **skewer** to remove Illumina adapters before aligning. Will run the shell script "skewer.sh".

```
#!/bin/bash
#SBATCH -D /mnt/lustre/macmaneslab/devon/bruce/fastq
#SBATCH -p macmanes,shared
#SBATCH --job-name="oro-skewer"
#SBATCH --ntasks=1
#SBATCH --output=skewer.log

module purge
module load linuxbrew/colsa

srun skewer -q 20 Bsuis0451_CAGAGAGG-TATCCTCT_L001_R1_001.fastq.gz Bsuis0451_CAGAGAGG-TATCCTCT_L001_R2_001.fastq.gz --output 0451-trimmed -z
```

