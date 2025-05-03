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
    """Create overlay plot of coadd and emode"""
    # Load and prepare both datasets
    coadd_hdul = fits.open(args.coaddpath)
    coadd_data = coadd_hdul[0].data
    coadd_header = coadd_hdul[0].header.copy()
    coadd_header['CTYPE1'] = 'RA---TAN'
    coadd_header['CTYPE2'] = 'DEC--TAN'
    wcs_coadd = WCS(coadd_header)
    
    emode_hdul = fits.open(args.emodepath)
    emode_data = emode_hdul[0].data
    emode_header = emode_hdul[0].header.copy()
    emode_header['CTYPE1'] = 'RA---TAN'
    emode_header['CTYPE2'] = 'DEC--TAN'
    wcs_emode = WCS(emode_header)

    # Crop both images around target coordinates
    def crop_image(data, wcs):
        target_x, target_y = wcs.world_to_pixel(target_coord)
        crop_range = 500
        height, width = data.shape
        y_start = max(0, int(target_y) - crop_range)
        y_end = min(height, int(target_y) + crop_range)
        x_start = max(0, int(target_x) - crop_range)
        x_end = min(width, int(target_x) + crop_range)
        return data[y_start:y_end, x_start:x_end]
    
    coadd_crop = crop_image(coadd_data, wcs_coadd)
    emode_crop = crop_image(emode_data, wcs_emode)

    # Create overlay plot using coadd's WCS
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs_coadd)
    
    # Plot coadd as background
    norm_coadd = ImageNormalize(coadd_crop, interval=ZScaleInterval())
    im = ax.imshow(coadd_crop, cmap='magma', norm=norm_coadd, origin='lower', alpha=0.8)

    # Plot emode as contours with auto-scaling
    norm_emode = ImageNormalize(emode_crop, interval=ZScaleInterval())
    levels = np.linspace(norm_emode.vmin, norm_emode.vmax, 7)
    contours = ax.contour(emode_crop, levels=levels, colors='cyan', linewidths=1.5, alpha=0.7)

    # Formatting
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

    # Add labels and legend
    plt.colorbar(im, pad=0.15).set_label('Coadd Flux (Jy/beam)')
    ax.clabel(contours, inline=True, fontsize=8, fmt='%1.2f')
    
    plt.savefig(os.path.join(args.outdir, f"{args.cluster_name}_{args.band_name}_overlay.png"))
    plt.close()
    coadd_hdul.close()
    emode_hdul.close()

if __name__ == '__main__':
    args = parse_args()
    target_coord = get_reference_coords(args.emodepath)
    
    # Generate individual plots
    process_image(args, args.coaddpath, 'coadd', target_coord)
    process_image(args, args.emodepath, 'emode', target_coord)
    
    # Generate overlay plot
    create_overlay(args, target_coord)