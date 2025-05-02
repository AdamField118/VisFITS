from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.wcs import WCS
from argparse import ArgumentParser
import os

def parse_args():

    parser = ArgumentParser()

    parser.add_argument('outdir', type=str, default=None,
                        help='Output directory for images')
    parser.add_argument('coaddpath', type=str, default=None,
                        help='Path to COADD File')
    parser.add_argument('emodepath', type=str, default=None,
                        help='Path to COADD File')
    parser.add_argument('cluster_name', type=str, default=None,
                        help='Path to COADD File')
    parser.add_argument('band_name', type=str, default=None,
                        help='Path to COADD File')

    return parser.parse_args()

def main(args, path, int):
    outdir = args.outdir
    cluster_name = args.cluster_name
    band_name = args.band_name

    hdul = fits.open(path)
    image_data = hdul[0].data
    header = hdul[0].header

    # Force TPV to be recognized as TAN projection (temporary workaround)
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'

    wcs = WCS(header)

    norm = ImageNormalize(image_data, interval=ZScaleInterval())

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs)

    center_x = int(header['CRPIX1'])
    center_y = int(header['CRPIX2'])
    range = 500

    # Plot a smaller section for testing
    im = ax.imshow(image_data[(center_y-range):(center_y+range), (center_x-range):(center_x+range)], cmap='magma', norm=norm, origin='lower')

    # Explicitly configure coordinates
    ra = ax.coords['ra']
    dec = ax.coords['dec']

    # International Celestial Reference System
    ra.set_axislabel('Right Ascension (ICRS)')
    dec.set_axislabel('Declination (ICRS)')

    # Force tick parameters
    ra.set_ticks(size=10, color='white', width=1.5)
    dec.set_ticks(size=10, color='white', width=1.5)

    # Format ticks in sexagesimal
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm:ss')

    # Add grid
    ax.coords.grid(True, color='white', linestyle='--', alpha=0.7)

    # Astronomical convention
    ax.invert_xaxis()


    # Janksy (Jy) unit of flux density: 1 Jy = 10^{-26} W/m^{2}/Hz
    # beam is the effective resolution element of the telescope
    # Jy/beam is the flux density per beam area (common unit in radio/mm astronomy)
    plt.colorbar(im, pad=0.15).set_label('Flux (Jy/beam)')
    plt.savefig(os.path.join(outdir, f"{cluster_name}_{band_name}_{int}.png"))
    hdul.close()

if __name__ == '__main__':
    args = parse_args()
    coaddpath = args.coaddpath
    emodepath = args.emodepath
    main(args, coaddpath, 1)
    main(args, emodepath, 2)