#!/bin/sh
#SBATCH -t 1:59:59
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=5g
#SBATCH --partition=short
#SBATCH -J VisFITS
#SBATCH -v
#SBATCH -o slurm_outfiles/out_visfits_%j.log
#SBATCH -e slurm_outfiles/err_visfits_%j.err

#Configuration

export COADDPATH="../superbit-lensing/data/PLCKG287d0p32d9/b/coadd/PLCKG287d0p32d9_coadd_b.fits"
export EMODEPATH="../SMPy/outputs/kaiser_squires/simulation_testing_kaiser_squires_e_mode.fits"
 
module load miniconda3
# Set Conda environment
export CONDA_ENV="visfits"

# Ensure the conda command is available
source ~/.bashrc
source activate $CONDA_ENV

dirname="slurm_outfiles"
if [ ! -d "$dirname" ]; then
     echo " Directory $dirname does not exist. Creating now"
     mkdir -p -- "$dirname"
     echo " $dirname created"
else
     echo " Directory $dirname exists"
fi

echo "Proceeding with code..."

python VisFITS.py "$OUTDIR" "$COADDPATH" "$EMODEPATH" "$cluster_name" "$band_name"

# Check if the Python script ran successfully
: '
if [ $? -eq 0 ]; then
    echo "Python script executed successfully. Check your outputs dir."
    exit 1
else
    echo "Python script failed. Check you slurm_outfiles dir."
    exit 1
fi
'