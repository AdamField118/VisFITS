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
    
    # Add this line to ensure proper WCS handling
    wcs = WCS(header).celestial  # Focus on celestial coordinates only

    # Create plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs)
    
    # Add these lines to explicitly set coordinate display
    ax.coords[0].set_axislabel('Right Ascension (ICRS)')
    ax.coords[1].set_axislabel('Declination (ICRS)')
    
    norm = ImageNormalize(image_data, interval=ZScaleInterval())
    im = ax.imshow(image_data, cmap='magma', norm=norm, origin='lower', aspect='auto')

    # Coordinate formatting - modified grid parameters
    ra = ax.coords[0]
    dec = ax.coords[1]
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm:ss')
    
    # Modified grid settings for better visibility
    ax.coords.grid(True, color='gray', linestyle='--', alpha=0.7)
    
    # Add this to ensure proper layout
    plt.tight_layout()
    
    plt.colorbar(im, pad=0.15).set_label('Flux (Jy/beam)')
    plt.savefig(os.path.join(args.outdir, f"{args.cluster_name}_{args.band_name}_{name}.png"), bbox_inches='tight')
    plt.close()
    hdul.close()

if __name__ == '__main__':
    args = parse_args()
    process_image(args, args.coaddpath, 'coadd')
    process_image(args, args.emodepath, 'emode')