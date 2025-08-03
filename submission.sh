#!/bin/sh
#SBATCH -t 1:59:59
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=8g
#SBATCH --partition=short
#SBATCH -J VisFITS_overlay
#SBATCH -v
#SBATCH -o slurm_outfiles/out_visfits_%j.log
#SBATCH -e slurm_outfiles/err_visfits_%j.err

# ============================================================================
# VisFITS Configuration - EDIT THESE PATHS FOR YOUR DATA
# ============================================================================

# REQUIRED: Full paths to your FITS files
COADDPATH="/home/adfield/weak_lensing/superbit-lensing/data/PLCKG287d0p32d9/b/coadd/PLCKG287d0p32d9_coadd_b.fits"
EMODEPATH="/home/adfield/weak_lensing/SMPy/outputs/kaiser_squires/simulation_testing_kaiser_squires_e_mode.fits"

# EXAMPLE paths (uncomment and modify for your setup):
# COADDPATH="../superbit-lensing/data/PLCKG287d0p32d9/b/coadd/PLCKG287d0p32d9_coadd_b.fits"
# EMODEPATH="../SMPy/outputs/kaiser_squires/simulation_testing_kaiser_squires_e_mode.fits"

# Optional: Customize output and display settings
OUTPUT_NAME="galaxy_cluster_analysis"  # Base name for output files
COADD_ALPHA="0.9"                     # Transparency for baryonic matter (0-1)
EMODE_ALPHA="0.6"                     # Transparency for dark matter (0-1)

# Conda environment name
CONDA_ENV="visfits"

# ============================================================================
# Script execution (usually don't need to modify below this line)
# ============================================================================

echo "=========================================="
echo "VisFITS: FITS Overlay Visualization"
echo "=========================================="

# Load conda environment
module load miniconda3
source ~/.bashrc
source activate $CONDA_ENV

# Create output directories
mkdir -p outputs slurm_outfiles

echo "Configuration:"
echo "  Coadd file: $COADDPATH"
echo "  E-mode file: $EMODEPATH"
echo "  Output name: $OUTPUT_NAME"
echo "  Coadd alpha: $COADD_ALPHA"
echo "  E-mode alpha: $EMODE_ALPHA"
echo ""

# Verify input files exist
if [ ! -f "$COADDPATH" ]; then
    echo "ERROR: Coadd file not found: $COADDPATH"
    echo "Please edit the COADDPATH variable in this script with the full path to your coadd FITS file"
    exit 1
fi

if [ ! -f "$EMODEPATH" ]; then
    echo "ERROR: E-mode file not found: $EMODEPATH"
    echo "Please edit the EMODEPATH variable in this script with the full path to your E-mode FITS file"
    exit 1
fi

echo "Input files verified ✓"
echo "Starting VisFITS processing..."

# Run VisFITS with all options
python visfits.py "$COADDPATH" "$EMODEPATH" \
    --output-name "$OUTPUT_NAME" \
    --coadd-alpha "$COADD_ALPHA" \
    --emode-alpha "$EMODE_ALPHA"

# Check result
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "SUCCESS!"
    echo "=========================================="
    echo ""
    echo "Generated files:"
    echo "  📊 ./outputs/${OUTPUT_NAME}_overlay.png   (Main overlay for poster)"
    echo "  📈 ./outputs/${OUTPUT_NAME}_process.png   (Process explanation)"
    echo ""
    echo "Main overlay shows:"
    echo "  🔘 Gray regions = Baryonic matter (visible galaxies)"
    echo "  🟣 Purple regions = Dark matter (from weak lensing)"
    echo ""
    echo "Process figure shows:"
    echo "  [Coadd] + [Mass Map] → [Overlay] visualization"
    echo ""
    echo "Perfect for your poster! The overlay reveals where"
    echo "dark matter is concentrated relative to visible matter."
else
    echo ""
    echo "❌ VisFITS processing failed!"
    echo "Check the error log: slurm_outfiles/err_visfits_*.err"
    exit 1
fi