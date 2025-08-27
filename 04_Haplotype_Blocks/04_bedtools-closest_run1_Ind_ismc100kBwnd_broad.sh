#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=80GB
#SBATCH --time=00:01:00
#SBATCH --account=ag-waldvogel
#SBATCH --job-name=bedtools_closest
#SBATCH --error /scratch/lpettric/jobs/%x-%N-%j.err
#SBATCH --output /scratch/lpettric/jobs/%x-%N-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lpettric@smail.uni-koeln.de

# Load modules
module purge
module load lang/Miniconda3/23.9.0-0
conda activate bedtools_env

OUTDIR=/projects/ag-waldvogel/pophistory/CRIP/comparison_HB_rho/run03

INPDIR=/projects/ag-waldvogel/pophistory/CRIP/plink_haplotype_blocks/run_250822

RHODIR=/projects/ag-waldvogel/pophistory/CRIP/ismc

IND=("MF1" "MF2" "MF3" "MF4" \
       "MG2" "MG3" "MG4" "MG5" \
       "NMF1" "NMF2" "NMF3" "NMF4" \
       "SI1" "SI2" "SI3" "SI4" \
       "SS1" "SS2" "SS3" "SS4")

CHR=("Chr1" "Chr2" "Chr3" "Chr4")

POP=("MF" "MG" "NMF" "SI" "SS")


# RHO: CONVERT BEDGRAPH TO BED AND SORT
### https://www.biostars.org/p/10141/
# Has been already done in 1.6

#for ind in "${IND[@]}"; do
#  for chr in "${CHR[@]}"; do    
#    cd ${RHODIR}/"${ind}"/"${chr}"
#    cat "${ind}"_"${chr}"_ismc.rho.100kb.bedgraph | grep -v '^chrom' \
#    | sed "s/chr1/${chr}/g" | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.100kb.bedgraph_renamed.bed
#    sort -k1,1 -k2,2n ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.100kb.bedgraph_renamed.bed > ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.100kb.bedgraph_renamed.sorted.bed
#  done
#done

wait

echo "Converting bedgraph to bed done"

wait

# NO MEAN RHO/10kb FOR EACH POP NEED If USING IndivAnalysis



# SPLIT HP-BED BY CHROMOSOME AND SORT
#for chr in "${CHR[@]}"; do
#    cd ${INPDIR}
#    awk 'NR>1 {printf "Chr%s\t%d\t%d\n", $1, $2-1, $3}' all_pop_blocks_broad.blocks.det > all_pop_haploblocks_broad.bed
#    cat all_pop_haploblocks_broad.bed | grep "${chr}" | sort -k1,1 -k2,2n > "${chr}"_all_pop_haploblocks_broad.bed
#done

#echo "Split HP-bed by chromosome and sort it done"

#wait

# BEDTOOLS CLOSEST
### with pophistory samples mapped on Crip4.0
### for each individual and each chromosome 
### HB compared to rho (1 kb)

cd ${OUTDIR}

for ind in "${IND[@]}"; do
  for chr in "${CHR[@]}"; do
    bedtools closest -d -t first -a ${RHODIR}/bed-files/"${ind}"_"${chr}"_ismc.rho.100kb.bedgraph_renamed.sorted.bed \
    -b ${INPDIR}/"${chr}"_all_pop_haploblocks_broad.bed > \
    "${ind}"_"${chr}"_ismc.rho.100kb.bedgraph_renamed.sorted_HB_broad.bed
  done
done

wait

echo "Bedtools closest finished"

wait

# FILTER FOR CLUSTERSIZE MORE OR EQUAL THAN 500 bp

for ind in "${IND[@]}"; do
  for chr in "${CHR[@]}"; do
    cd ${OUTDIR}
    awk '$7 - $6 >= 500' "${ind}"_"${chr}"_ismc.rho.100kb.bedgraph_renamed.sorted_HB_broad.bed > "${ind}"_"${chr}"_ismc.rho.100kb.bedgraph_renamed.sorted_HB_more500bp_broad.bed
  done
done

# FILTER FOR CLUSTERSIZE LESS THAN 500 bp

for ind in "${IND[@]}"; do
  for chr in "${CHR[@]}"; do
    cd ${OUTDIR}
    awk '$7 - $6 < 500' "${ind}"_"${chr}"_ismc.rho.100kb.bedgraph_renamed.sorted_HB_broad.bed > "${ind}"_"${chr}"_ismc.rho.100kb.bedgraph_renamed.sorted_HB_less500bp_broad.bed
  done
done

