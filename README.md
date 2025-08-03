# VisFITS: Dark Matter vs Baryonic Matter Overlay

Creates overlay visualizations of weak lensing data to show dark matter distribution relative to baryonic matter around galaxy clusters.

## Quick Start

1. **Set up environment:**
   ```bash
   conda env create -f visfits.yml
   ```

2. **Edit file paths in `submission.sh`:**
   ```bash
   COADDPATH="/path/to/your/coadd.fits"
   EMODEPATH="/path/to/your/emode.fits"
   ```

3. **Run:**
   ```bash
   sbatch submission.sh
   ```

4. **Get your overlay:**
   ```bash
   ls outputs/
   # -> galaxy_cluster_analysis_overlay.png
   ```

## What You Get

**One overlay image** showing:
- **Gray regions:** Baryonic matter (visible galaxies from coadd)
- **Purple regions:** Dark matter distribution (from E-mode reconstruction)
**One flow chart** showing:
- Coadd image
- E-mode mass map
- Then the same overlay image

## Requirements

- Completed SuperBIT lensing pipeline (coadd FITS file)
- SMPy E-mode reconstruction (E-mode FITS file)
- SLURM cluster access

## Customization

Edit these variables in `submission.sh`:
- `OUTPUT_NAME` - Base name for output files  
- `COADD_ALPHA` - Transparency of baryonic matter layer (0-1)
- `EMODE_ALPHA` - Transparency of dark matter layer (0-1)

## Output

Both files in the `./outputs/` directory.

This overlay reveals where dark matter is concentrated relative to visible matter.