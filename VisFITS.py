from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, ZScaleInterval
from astropy.wcs import WCS
from argparse import ArgumentParser
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

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
    
    # Force consistent projection
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    wcs = WCS(header)
    
    # Convert target coordinates to this image's pixel space
    target_x, target_y = wcs.world_to_pixel(target_coord)
    
    crop_range = 500
    height, width = image_data.shape
    
    # Calculate safe crop boundaries
    y_start = max(0, int(target_y) - crop_range)
    y_end = min(height, int(target_y) + crop_range)
    x_start = max(0, int(target_x) - crop_range)
    x_end = min(width, int(target_x) + crop_range)
    
    cropped_data = image_data[y_start:y_end, x_start:x_end]

    # Create plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs)
    norm = ImageNormalize(cropped_data, interval=ZScaleInterval())
    im = ax.imshow(cropped_data, cmap='magma', norm=norm, origin='lower')

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
    ax.invert_xaxis()

    plt.colorbar(im, pad=0.15).set_label('Flux (Jy/beam)')
    plt.savefig(os.path.join(args.outdir, f"{args.cluster_name}_{args.band_name}_{name}.png"))
    plt.close()
    hdul.close()

def create_overlay(args, target_coord):
    """Create overlay plot with properly scaled contours"""
    # Load data with context managers
    with fits.open(args.coaddpath) as coadd_hdul, fits.open(args.emodepath) as emode_hdul:
        # Process coadd data
        coadd_data = coadd_hdul[0].data
        coadd_header = coadd_hdul[0].header.copy()
        coadd_header['CTYPE1'] = 'RA---TAN'
        coadd_header['CTYPE2'] = 'DEC--TAN'
        wcs_coadd = WCS(coadd_header)

        # Process emode data
        emode_data = emode_hdul[0].data
        emode_header = emode_hdul[0].header.copy()
        emode_header['CTYPE1'] = 'RA---TAN'
        emode_header['CTYPE2'] = 'DEC--TAN'
        wcs_emode = WCS(emode_header)

        # Crop both images using common WCS
        def get_cropped(data, wcs):
            target_x, target_y = wcs.world_to_pixel(target_coord)
            crop_range = 500
            y_start = max(0, int(target_y) - crop_range)
            y_end = min(data.shape[0], int(target_y) + crop_range)
            x_start = max(0, int(target_x) - crop_range)
            x_end = min(data.shape[1], int(target_x) + crop_range)
            return data[y_start:y_end, x_start:x_end]

        coadd_crop = get_cropped(coadd_data, wcs_coadd)
        emode_crop = get_cropped(emode_data, wcs_emode)

        # Create figure with coadd WCS
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection=wcs_coadd)

        # Plot coadd with adjusted transparency
        norm_coadd = ImageNormalize(coadd_crop, interval=ZScaleInterval())
        im = ax.imshow(coadd_crop, cmap='magma', norm=norm_coadd, origin='lower', alpha=0.7)

        # Calculate smart contour levels for emode
        emode_abs_max = np.nanmax(np.abs(emode_crop))
        contour_levels = np.linspace(-emode_abs_max, emode_abs_max, 11)

        # Plot contours with enhanced visibility
        contours = ax.contour(
            emode_crop, 
            levels=contour_levels,
            colors=['limegreen' if v > 0 else 'magenta' for v in contour_levels],
            linewidths=1.2,
            linestyles='-',
            alpha=0.9
        )

        # Add contour labels
        ax.clabel(contours, inline=True, fontsize=8, fmt='%1.2e', colors='white')

        # Add legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], color='limegreen', lw=2, label='Positive E-mode'),
            Line2D([0], [0], color='magenta', lw=2, label='Negative E-mode')
        ]
        ax.legend(handles=legend_elements, loc='upper right', 
                facecolor='black', edgecolor='white',
                fontsize=9, labelcolor='white')

        # Format axes (same as before)
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

        # Add colorbar
        cbar = plt.colorbar(im, pad=0.15)
        cbar.set_label('Coadd Flux (Jy/beam)')

        plt.savefig(os.path.join(args.outdir, f"{args.cluster_name}_{args.band_name}_overlay.png"))
        plt.close()

if __name__ == '__main__':
    args = parse_args()
    target_coord = get_reference_coords(args.emodepath)
    
    # Generate individual plots
    process_image(args, args.coaddpath, 'coadd', target_coord)
    process_image(args, args.emodepath, 'emode', target_coord)
    
    # Generate overlay plot
    create_overlay(args, target_coord)