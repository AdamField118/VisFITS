#!/bin/sh
#SBATCH -t 1:59:59
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=5g
#SBATCH --partition=short
#SBATCH -J VizFITS
#SBATCH -v
#SBATCH -o slurm_outfiles/out_vizfits_%j.log
#SBATCH -e slurm_outfiles/err_vizfits_%j.err

# Load configuration file
source "$SLURM_SUBMIT_DIR/config.sh"

dirname="slurm_outfiles"
if [ ! -d "$dirname" ]; then
     echo " Directory $dirname does not exist. Creating now"
     mkdir -p -- "$dirname"
     echo " $dirname created"
else
     echo " Directory $dirname exists"
fi

echo "Proceeding with code..."

python VisFITS.py "$OUTDIR" "$COADDPATH" "$cluster_name" "$band_name"

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