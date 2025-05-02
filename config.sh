# config.sh - Configuration file for submission.sh

export cluster_name="PLCKG287d0p32d9"
export band_name="b"

export DATADIR="./data/"
export OUTDIR="./outputs/"
export COADDPATH="${DATADIR}/${cluster_name}/coadd/${band_name}.fits"
 
module load miniconda3
# Set Conda environment
export CONDA_ENV="visfits"

# Ensure the conda command is available
source ~/.bashrc
source activate $CONDA_ENV