#!/bin/bash
#SBATCH --job-name=nf_test    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=benjamin.clark.hsj@ssss.gouv.qc.ca     # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=2G                     # Job memory request
#SBATCH --time=05:00:00               # Time limit hrs:min:sec
#SBATCH --output=nf_test_%j.log   # Standard output and error log
#SBATCH --account=rrg-jacquese        #group account

module load nextflow

nextflow run main.nf -c conf/ccdb.config nextflow.config -with-report -resume

