#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=5
#SBATCH --mem=150GB
#SBATCH --time=200:00:00
#SBATCH --account=ag-waldvogel
#SBATCH --job-name=crip_plink-r2-blocks
#SBATCH --error /scratch/lpettric/jobs/%x-%N-%j.err
#SBATCH --output /scratch/lpettric/jobs/%x-%N-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lpettric@smail.uni-koeln.de

# Define paths and variables
IND=("MF1" "MF2" "MF3" "MF4" \
       "MG2" "MG3" "MG4" "MG5" \
       "NMF1" "NMF2" "NMF3" "NMF4" \
       "SI1" "SI2" "SI3" "SI4" \
       "SS1" "SS2" "SS3" "SS4")

CHR=("Chr1" "Chr2" "Chr3" "Chr4")

POP=("MF" "MG" "NMF" "SI" "SS")

MF=("MF1" "MF2" "MF3" "MF4")

MG=("MG2" "MG3" "MG4" "MG5")

NMF=("NMF1" "NMF2" "NMF3" "NMF4")

SI=("SI1" "SI2" "SI3" "SI4")

SS=("SS1" "SS2" "SS3" "SS4")

OUTDIR=/projects/ag-waldvogel/pophistory/CRIP/plink_haplotype_blocks/run_250822


# Load modules
module purge
module load bio/BCFtools/1.19-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load tools/parallel/20240322-GCCcore-13.2.0
module load lang/Miniconda3/23.9.0-0
conda activate plink_env

# OpenAI ChatGPT (version 4), assistance in coding. Available at: https://openai.com/chatgpt

# Merge VCF files per population
# Convert merged VCF files to PLINK format for each population
cd $OUTDIR

parallel -j 5 plink --bfile {}_filtered   --blocks no-pheno-req no-small-max-span --blocks-max-kb 1000 --blocks-min-maf 0.10 --blocks-strong-lowci 0.65 --blocks-strong-highci 0.97 --blocks-inform-frac 0.90 --snps-only --threads 5 --out {}_blocks ::: "${POP[@]}"

parallel -j 2 plink --bfile {}_filtered --r2 --threads 5 --out {}_filtered_ld ::: "${POP[@]}"

echo "PLINK analysis for each population completed"

# NOT WORKING
# Perform LD Calculation and Haplotype Block Analysis for all populations combined
plink --bfile all_pop_filtered --blocks no-pheno-req no-small-max-span --blocks-max-kb 1000 --blocks-min-maf 0.10 --blocks-strong-lowci 0.65 --blocks-strong-highci 0.97 --blocks-inform-frac 0.90 --snps-only --threads 5 --out all_pop_blocks_broad
  
plink --bfile all_pop_filtered --r2 --threads 5 --out all_pop_filtered_ld

echo "PLINK analysis for combined populations completed"
echo "PLINK analysis completed"


# read .blocks.det in R
