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
    
    wcs = WCS(header).celestial
    print(wcs.to_header())

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs)
    
    lon = ax.coords[0]
    lat = ax.coords[1]
    
    lon.set_axislabel('Right Ascension (ICRS)', minpad=0.8)
    lon.set_major_formatter('hh:mm:ss')
    lon.set_ticks(color='black', size=10, width=1.5)
    lon.display_minor_ticks(True)
    
    lat.set_axislabel('Declination (ICRS)', minpad=0.9)
    lat.set_major_formatter('dd:mm:ss')
    lat.set_ticks(color='black', size=10, width=1.5)
    lat.display_minor_ticks(True)
    
    ax.set_xlim(-0.5, header['NAXIS1']-0.5)
    ax.set_ylim(-0.5, header['NAXIS2']-0.5)
    
    norm = ImageNormalize(image_data, interval=ZScaleInterval())
    im = ax.imshow(image_data, cmap='magma', norm=norm, 
                  origin='lower', aspect='auto')

    ax.coords.grid(True, color='gray', linestyle='--', alpha=0.7)
    
    cbar = plt.colorbar(im, pad=0.05)
    cbar.set_label('Flux (Jy/beam)', rotation=270, labelpad=25)

    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)
    plt.savefig(os.path.join(args.outdir, f"{args.cluster_name}_{args.band_name}_{name}.png"), 
                bbox_inches='tight', dpi=150)
    plt.close()
    hdul.close()

if __name__ == '__main__':
    args = parse_args()
    process_image(args, args.coaddpath, 'coadd')
    process_image(args, args.emodepath, 'emode')