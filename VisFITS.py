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
    
    # Handle WCS with distortion-aware parsing
    wcs = WCS(header, relax=True).celestial

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs)

    # Get pixel scales using Astropy's built-in method
    try:
        pixel_scales = wcs.proj_plane_pixel_scales() * 3600  # Convert to arcsec
        cdelt1, cdelt2 = abs(pixel_scales[0].value), abs(pixel_scales[1].value)
    except AttributeError:
        # Fallback for older Astropy versions
        cdelt1 = abs(wcs.pixel_scale_matrix[0,0] * 3600)
        cdelt2 = abs(wcs.pixel_scale_matrix[1,1] * 3600)

    # Calculate spans using header dimensions
    ra_span = cdelt1 * header['NAXIS1'] 
    dec_span = cdelt2 * header['NAXIS2']

    # Dynamic formatting based on actual span
    ax.coords[0].set_axislabel('Right Ascension (ICRS)')
    ax.coords[1].set_axislabel('Declination (ICRS)')
    
    # Set formatters based on angular size
    if ra_span < 60:  # Less than 1 arcminute
        ax.coords[0].set_major_formatter('hh:mm:ss.ss')
    else:
        ax.coords[0].set_major_formatter('hh:mm:ss')

    if dec_span < 60:
        ax.coords[1].set_major_formatter('dd:mm:ss.ss')
    else:
        ax.coords[1].set_major_formatter('dd:mm:ss')

    # Handle TPV-specific issues
    if 'TPV' in header.get('CTYPE1', ''):
        ax.coords[0].display_minor_ticks(False)
        ax.coords[1].display_minor_ticks(False)
        grid_style = dict(color='gray', linestyle=':', alpha=0.5)
    else:
        grid_style = dict(color='gray', linestyle='--', alpha=0.7)
    
    ax.coords.grid(True, **grid_style)

    # Set axis limits
    ax.set_xlim(-0.5, header['NAXIS1'] - 0.5)
    ax.set_ylim(-0.5, header['NAXIS2'] - 0.5)

    # Plot image
    norm = ImageNormalize(image_data, interval=ZScaleInterval())
    im = ax.imshow(image_data, cmap='magma', norm=norm,
                  origin='lower', aspect='auto', interpolation='nearest')

    # Add colorbar
    cbar = plt.colorbar(im, pad=0.05)
    cbar.set_label('Flux (Jy/beam)', rotation=270, labelpad=25)

    # Save output
    plt.savefig(os.path.join(args.outdir, f"{args.cluster_name}_{args.band_name}_{name}.png"),
                bbox_inches='tight', dpi=150)
    plt.close()
    hdul.close()

if __name__ == '__main__':
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    process_image(args, args.coaddpath, 'coadd')
    process_image(args, args.emodepath, 'emode')