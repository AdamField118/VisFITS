from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.wcs import WCS
from argparse import ArgumentParser
import os
from astropy.coordinates import SkyCoord
import astropy.units as u

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('outdir', type=str, help='Output directory for images')
    parser.add_argument('coaddpath', type=str, help='Path to COADD File')
    parser.add_argument('emodepath', type=str, help='Path to EMODE File')
    parser.add_argument('cluster_name', type=str, help='Cluster name')
    parser.add_argument('band_name', type=str, help='Band name')
    return parser.parse_args()

def get_reference_coords(emode_path):
    """Get center coordinates from emode FITS header"""
    with fits.open(emode_path) as hdul:
        header = hdul[0].header
        return SkyCoord(header['CRVAL1'] * u.deg, 
                       header['CRVAL2'] * u.deg, 
                       frame='icrs')

def process_image(args, path, name, target_coord):
    hdul = fits.open(path)
    image_data = hdul[0].data
    header = hdul[0].header.copy()
    print(header)
    
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    wcs = WCS(header)

    # Create plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs)
    norm = ImageNormalize(image_data, interval=ZScaleInterval())
    im = ax.imshow(image_data, cmap='magma', norm=norm, origin='lower')

    # Coordinate formatting
    ra = ax.coords['ra']
    dec = ax.coords['dec']
    ra.set_axislabel('Right Ascension (ICRS)')
    dec.set_axislabel('Declination (ICRS)')
    ra.set_ticks(size=10, color='white', width=1.5)
    dec.set_ticks(size=10, color='white', width=1.5)
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm:ss')
    ax.coords.grid(True, color='white', linestyle='--', alpha=0.7)

    plt.colorbar(im, pad=0.15).set_label('Flux (Jy/beam)')
    plt.savefig(os.path.join(args.outdir, f"{args.cluster_name}_{args.band_name}_{name}.png"))
    plt.close()
    hdul.close()

if __name__ == '__main__':
    args = parse_args()
    
    # Get reference coordinates from emode
    target_coord = get_reference_coords(args.emodepath)
    
    # Process both images using the same reference coordinates
    process_image(args, args.coaddpath, 'coadd', target_coord)
    process_image(args, args.emodepath, 'emode', target_coord)