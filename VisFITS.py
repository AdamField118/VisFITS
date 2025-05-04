from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.wcs import WCS
from argparse import ArgumentParser
import os

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('outdir', type=str, help='Output directory for images')
    parser.add_argument('coaddpath', type=str, help='Path to COADD File')
    parser.add_argument('emodepath', type=str, help='Path to EMODE File')
    parser.add_argument('cluster_name', type=str, help='Cluster name')
    parser.add_argument('band_name', type=str, help='Band name')
    return parser.parse_args()

def process_image(args, path, name):
    hdul = fits.open(path)
    image_data = hdul[0].data
    header = hdul[0].header
    
    # Handle TPV distortion explicitly
    wcs = WCS(header, relax=True).celestial  # Force celestial axes

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs)
    
    # Universal settings for both projections
    ax.coords[0].set_axislabel('Right Ascension (ICRS)', minpad=0.8)
    ax.coords[1].set_axislabel('Declination (ICRS)', minpad=0.8)
    
    # Dynamic formatting based on actual coordinate span
    ra_span = abs(header['CDELT1'] * header['NAXIS1'] * 3600)  # arcsec
    dec_span = abs(header['CDELT2'] * header['NAXIS2'] * 3600)
    
    # Auto-adjust tick density
    if ra_span < 60:  # Small RA span (< 1 arcmin)
        ax.coords[0].set_major_formatter('hh:mm:ss.s')
    else:
        ax.coords[0].set_major_formatter('hh:mm:ss')
        
    if dec_span < 60:  # Small Dec span
        ax.coords[1].set_major_formatter('dd:mm:ss.s')
    else:
        ax.coords[1].set_major_formatter('dd:mm:ss')
    
    # TPV-specific fixes
    if 'TPV' in header['CTYPE1']:
        ax.coords[0].display_minor_ticks(False)
        ax.coords[1].display_minor_ticks(False)
        ax.coords.grid(True, color='gray', linestyle=':', alpha=0.5)
    else:
        ax.coords.grid(True, color='gray', linestyle='--', alpha=0.7)

    # Dynamic axis scaling
    ax.set_xlim(-0.5, header['NAXIS1'] - 0.5)
    ax.set_ylim(-0.5, header['NAXIS2'] - 0.5)
    
    # Image plotting
    norm = ImageNormalize(image_data, interval=ZScaleInterval())
    im = ax.imshow(image_data, cmap='magma', norm=norm, 
                  origin='lower', aspect='auto', interpolation='nearest')

    # Colorbar with dynamic padding
    cbar = plt.colorbar(im, pad=0.05)
    cbar.set_label('Flux (Jy/beam)', rotation=270, labelpad=25)

    plt.savefig(os.path.join(args.outdir, f"{args.cluster_name}_{args.band_name}_{name}.png"),
                bbox_inches='tight', dpi=150)
    plt.close()
    hdul.close()

if __name__ == '__main__':
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    process_image(args, args.coaddpath, 'coadd')
    process_image(args, args.emodepath, 'emode')