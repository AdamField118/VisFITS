# config.sh - Configuration file for submission.sh

export cluster_name="PLCKG287d0p32d9"
export band_name="b"
export method="kaiser_squires"

export DATADIR="../superbit-lensing/data"
export OUTDIR="./outputs"
export COADDPATH="${DATADIR}/${cluster_name}/${band_name}/coadd/${cluster_name}_coadd_${band_name}.fits"
export EMODEPATH="../SMPy/outputs/${method}/simulation_testing_${method}_e_mode.fits"
 
module load miniconda3
# Set Conda environment
export CONDA_ENV="visfits"

# Ensure the conda command is available
source ~/.bashrc
source activate $CONDA_ENV