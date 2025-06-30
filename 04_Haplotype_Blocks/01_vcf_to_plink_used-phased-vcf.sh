#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=150GB
#SBATCH --time=72:00:00
#SBATCH --account=ag-waldvogel
#SBATCH --job-name=crip_plink-r2-blocks
#SBATCH --error /scratch/lpettric/jobs/%x-%N-%j.err
#SBATCH --output /scratch/lpettric/jobs/%x-%N-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lpettric@smail.uni-koeln.de

# Load required modules
module purge
module load bio/BCFtools/1.19-GCC-13.2.0
module load bio/SAMtools/1.19.2-GCC-13.2.0
module load lang/Miniconda3/23.9.0-0
conda activate plink_env

# Define paths and variables
IND=("MF1" "MF2" "MF3" "MF4"
     "MG2" "MG3" "MG4" "MG5"
     "NMF1" "NMF2" "NMF3" "NMF4"
     "SI1" "SI2" "SI3" "SI4"
     "SS1" "SS2" "SS3" "SS4")

POP=("MF" "MG" "NMF" "SI" "SS")

INPUT_BASE="/projects/ag-waldvogel/pophistory/CRIP/phasing"
OUTDIR="/projects/ag-waldvogel/pophistory/CRIP/plink_haplotype_blocks"
FILEEND="_phased_merged.vcf.gz"
LOGFILE="${OUTDIR}/pipeline.log"

mkdir -p "$OUTDIR"
exec > >(tee -a "$LOGFILE") 2>&1  # Log stdout & stderr

echo "Pipeline started at $(date)"

# Step 1: Merge chromosomes per individual
echo "Merging chromosomes for each individual..."
for i in "${IND[@]}"; do
    INPUT_FILES=()
    for chr in {1..4}; do
        VCF_FILE="${INPUT_BASE}/Chr${chr}/phased/${i}_Chr${chr}${FILEEND}"
        if [[ -f "$VCF_FILE" ]]; then
            echo "Indexing $VCF_FILE"
            tabix -f -p vcf "$VCF_FILE"  # Ensure input VCF files are indexed
            
            # Assign unique variant IDs before merging
            FIXED_VCF="${OUTDIR}/${i}_Chr${chr}_fixed.vcf.gz"
            bcftools annotate --set-id '%CHROM:%POS' -Oz -o "$FIXED_VCF" "$VCF_FILE"
            tabix -p vcf "$FIXED_VCF"
            
            INPUT_FILES+=("$FIXED_VCF")
        else
            echo "WARNING: Missing file $VCF_FILE, skipping..."
        fi
    done

    if [[ ${#INPUT_FILES[@]} -gt 0 ]]; then
        echo "Merging ${INPUT_FILES[@]} for ${i}"
        bcftools concat "${INPUT_FILES[@]}" -Oz -o "${OUTDIR}/${i}_allChr.vcf.gz"
        
        # Reheader to correct the sample names (remove prefixes)
        bcftools reheader -s <(echo "$i") "${OUTDIR}/${i}_allChr.vcf.gz" -o "${OUTDIR}/${i}_allChr_cleaned.vcf.gz"
        
        # Index the reheadered VCF
        tabix -p vcf "${OUTDIR}/${i}_allChr_cleaned.vcf.gz"
        
        echo "Merged and cleaned VCF for ${i}: ${OUTDIR}/${i}_allChr_cleaned.vcf.gz"
    else
        echo "ERROR: No VCF files found for $i. Skipping..."
    fi
done

# Step 2: Merge individuals within each population
echo "Merging individuals within each population..."
for p in "${POP[@]}"; do
    POP_IND=()
    for i in "${IND[@]}"; do
        if [[ $i == $p* ]]; then
            POP_IND+=("${OUTDIR}/${i}_allChr_cleaned.vcf.gz")
        fi
    done

    if [[ ${#POP_IND[@]} -gt 0 ]]; then
        echo "Merging ${POP_IND[@]} for population ${p}"
        bcftools merge "${POP_IND[@]}" -Oz -o "${OUTDIR}/${p}_merged.vcf.gz"
        tabix -p vcf "${OUTDIR}/${p}_merged.vcf.gz"

        # Assign unique variant IDs to population-level VCF
        bcftools annotate --set-id '%CHROM:%POS' -Oz -o "${OUTDIR}/${p}_merged_fixed.vcf.gz" "${OUTDIR}/${p}_merged.vcf.gz"
        tabix -p vcf "${OUTDIR}/${p}_merged_fixed.vcf.gz"
        mv "${OUTDIR}/${p}_merged_fixed.vcf.gz" "${OUTDIR}/${p}_merged.vcf.gz"  # Replace old file
    else
        echo "ERROR: No merged individual VCFs found for $p. Skipping..."
    fi
done

# Step 3: Filter for SNPs only
echo "Filtering for SNPs only..."
for p in "${POP[@]}"; do
    if [[ -f "${OUTDIR}/${p}_merged.vcf.gz" ]]; then
        bcftools view -v snps -Oz -o "${OUTDIR}/${p}_snps.vcf.gz" "${OUTDIR}/${p}_merged.vcf.gz"
        tabix -p vcf "${OUTDIR}/${p}_snps.vcf.gz"
    else
        echo "WARNING: ${OUTDIR}/${p}_merged.vcf.gz not found, skipping SNP filtering..."
    fi
done

# Step 4: Convert to PLINK format
echo "Converting to PLINK format..."
for p in "${POP[@]}"; do
    if [[ -f "${OUTDIR}/${p}_snps.vcf.gz" ]]; then
        plink --vcf "${OUTDIR}/${p}_snps.vcf.gz" --make-bed --out "${OUTDIR}/${p}_plink" --threads 5
    else
        echo "WARNING: ${OUTDIR}/${p}_snps.vcf.gz not found, skipping PLINK conversion..."
    fi
done

# Step 5: Merge all populations into a single dataset
# Step 5: Merge all population VCF files into a single dataset
echo "Merging all population VCF files into a single dataset..."

# Create an array to store existing VCF file paths
VCF_FILES=()

# Collect all valid VCF files
for p in "${POP[@]}"; do
    VCF_PATH="${OUTDIR}/${p}_snps.vcf.gz"
    if [[ -f "$VCF_PATH" ]]; then
        VCF_FILES+=("$VCF_PATH")
    else
        echo "WARNING: $VCF_PATH not found, skipping..."
    fi
done

# Check if we have at least two files to merge
if [[ ${#VCF_FILES[@]} -ge 2 ]]; then
    # Merge all valid VCFs in one go
    bcftools merge "${VCF_FILES[@]}" -Oz -o "${OUTDIR}/all_pop_snps.vcf.gz"

    # Index the merged VCF
    tabix -p vcf "${OUTDIR}/all_pop_snps.vcf.gz"

    echo "VCF merging completed: ${OUTDIR}/all_pop_snps.vcf.gz"
else
    echo "ERROR: Not enough valid VCF files to merge. Merging skipped."
fi


# Convert merged VCF to PLINK format
echo "Converting merged VCF to PLINK format..."
plink --vcf "${OUTDIR}/all_pop_snps.vcf.gz" --make-bed --out "${OUTDIR}/all_pop_plink"


echo "Pipeline completed successfully at $(date)"
#end
