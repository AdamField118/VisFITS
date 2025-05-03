from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.wcs import WCS
from argparse import ArgumentParser
import os
import numpy as np
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

def load_and_prepare(path):
    """New helper function to load FITS data"""
    hdul = fits.open(path)
    data = hdul[0].data
    header = hdul[0].header.copy()
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    return data, header, hdul

def main(args, path, name):
    data, header, hdul = load_and_prepare(path)  # Modified line
    wcs = WCS(header)
    center_x = int(header['CRPIX1']) - 1
    center_y = int(header['CRPIX2']) - 1
    crop_range = 500
    cropped_data = data[
        (center_y - crop_range):(center_y + crop_range),
        (center_x - crop_range):(center_x + crop_range)
    ]

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs)
    norm = ImageNormalize(cropped_data, interval=ZScaleInterval())
    im = ax.imshow(cropped_data, cmap='magma', norm=norm, origin='lower')

    # Rest of your original plotting code remains unchanged
    ra = ax.coords['ra']
    dec = ax.coords['dec']
    ra.set_axislabel('Right Ascension (ICRS)')
    dec.set_axislabel('Declination (ICRS)')
    ra.set_ticks(size=10, color='white', width=1.5)
    dec.set_ticks(size=10, color='white', width=1.5)
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm:ss')
    ax.coords.grid(True, color='white', linestyle='--', alpha=0.7)
    ax.invert_xaxis()

    plt.colorbar(im, pad=0.15).set_label('Flux (Jy/beam)')
    plt.savefig(os.path.join(args.outdir, f"{args.cluster_name}_{args.band_name}_{name}.png"))
    plt.close()
    hdul.close()

def create_overlay(args):
    """Create overlay plot with adjusted scaling"""
    coadd_data, coadd_header, coadd_hdul = load_and_prepare(args.coaddpath)
    emode_data, emode_header, emode_hdul = load_and_prepare(args.emodepath)
    
    # Get sky center from coadd
    center_sky = SkyCoord(coadd_header['CRVAL1'] * u.deg, 
                        coadd_header['CRVAL2'] * u.deg, frame='icrs')
    
    # Convert to pixel coordinates in both images
    wcs_coadd = WCS(coadd_header)
    wcs_emode = WCS(emode_header)
    x_coadd, y_coadd = wcs_coadd.world_to_pixel(center_sky)
    x_emode, y_emode = wcs_emode.world_to_pixel(center_sky)

    # Crop both images
    crop_range = 500
    coadd_crop = coadd_data[
        int(y_coadd)-crop_range : int(y_coadd)+crop_range,
        int(x_coadd)-crop_range : int(x_coadd)+crop_range
    ]
    emode_crop = emode_data[
        int(y_emode)-crop_range : int(y_emode)+crop_range,
        int(x_emode)-crop_range : int(x_emode)+crop_range
    ]

    # Create overlay plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs_coadd)
    
    # Plot coadd with its own scale
    norm_coadd = ImageNormalize(coadd_crop, interval=ZScaleInterval())
    im = ax.imshow(coadd_crop, cmap='magma', norm=norm_coadd, origin='lower', alpha=0.7)

    # Calculate smart contours for emode
    emode_abs_max = np.nanmax(np.abs(emode_crop))
    positive_levels = np.linspace(0.1*emode_abs_max, emode_abs_max, 5)
    negative_levels = np.linspace(-0.1*emode_abs_max, -emode_abs_max, 5)
    
    # Plot positive and negative contours separately with different colors
    positive_contours = ax.contour(emode_crop, levels=positive_levels, 
                                 colors='limegreen', linewidths=1.5, linestyles='solid')
    negative_contours = ax.contour(emode_crop, levels=negative_levels, 
                                 colors='magenta', linewidths=1.5, linestyles='dashed')

    # Add contour labels
    ax.clabel(positive_contours, inline=True, fontsize=8, fmt='%1.2f')
    ax.clabel(negative_contours, inline=True, fontsize=8, fmt='%1.2f')

    # Formatting
    cbar = plt.colorbar(im, pad=0.15)
    cbar.set_label('Coadd Flux (Jy/beam)')
    
    # Add contour legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='limegreen', lw=2, label='Positive E-mode'),
        Line2D([0], [0], color='magenta', lw=2, linestyle='--', label='Negative E-mode')
    ]
    ax.legend(handles=legend_elements, loc='upper right', facecolor='black', 
            edgecolor='white', fontsize=10, labelcolor='white')

    # Keep existing coordinate formatting
    ra = ax.coords['ra']
    dec = ax.coords['dec']
    ra.set_axislabel('Right Ascension (ICRS)')
    dec.set_axislabel('Declination (ICRS)')
    ra.set_ticks(size=10, color='white', width=1.5)
    dec.set_ticks(size=10, color='white', width=1.5)
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm:ss')
    ax.coords.grid(True, color='white', linestyle='--', alpha=0.7)
    ax.invert_xaxis()

    plt.savefig(os.path.join(args.outdir, f"{args.cluster_name}_{args.band_name}_overlay.png"))
    plt.close()
    coadd_hdul.close()
    emode_hdul.close()

if __name__ == '__main__':
    args = parse_args()
    main(args, args.coaddpath, 'coadd')  # Original call
    main(args, args.emodepath, 'emode')  # Original call
    create_overlay(args)  # New call for overlay