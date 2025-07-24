#!/bin/bash
#SBATCH --job-name=phylostrat
#SBATCH --array=33-41
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=3-00:00
#SBATCH --output=/scratch/sgpl222/stratigraphy/logs/phylostrat_%A_%a.out
#SBATCH --error=/scratch/sgpl222/stratigraphy/logs/phylostrat_%A_%a.err
#SBATCH --account=coa_bmjo263_uksr


module purge
source ~/.bashrc
conda activate blast_env
module load ccs/java/jdk1.8.0_202

INPUT_DIR="/scratch/sgpl222/stratigraphy/formatted_proteins"
OUTPUT_DIR="/scratch/sgpl222/stratigraphy/phylostrat_outputs"
SCRIPT="/scratch/sgpl222/stratigraphy/createPSmap.pl"
DB="/scratch/sgpl222/stratigraphy/protein_db.faa"
SPECIES_LIST="/scratch/sgpl222/stratigraphy/species_list.txt"

species=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SPECIES_LIST")

echo "Running createPSmap.pl for $species"
echo "Input file: $INPUT_DIR/${species}.faa"
echo "Database: $DB"
echo "Output prefix: $OUTPUT_DIR/${species}_animals"

mkdir -p "$OUTPUT_DIR"

$SCRIPT \
    -o "$INPUT_DIR/${species}.faa" \
    --database "$DB" \
    --prefix "${species}_animals" \
    --seqOffset 50 \
    --evalue 1e-5 \
    --threads 12  \
    --blastPlus


