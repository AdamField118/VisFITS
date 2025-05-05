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
    parser.add_argument('coaddpath', type=str, help='Path to COADD File')
    parser.add_argument('emodepath', type=str, help='Path to EMODE File')
    return parser.parse_args()

def get_reference_coords(emode_path):
    """Get center coordinates from emode FITS header"""
    with fits.open(emode_path) as hdul:
        header = hdul[0].header
        return SkyCoord(header['CRVAL1'] * u.deg,
                       header['CRVAL2'] * u.deg,
                       frame='icrs')

def configure_wcs(original_header, target_coord):
    """Create aligned WCS while preserving original projection"""
    # Create a copy of the header to modify
    header = original_header.copy()
    
    # Update reference coordinates while keeping original projection
    header['CRVAL1'] = target_coord.ra.deg
    header['CRVAL2'] = target_coord.dec.deg
    
    # Handle both SIP (TPV) and TAN conventions
    wcs = WCS(header, relax=True)
    
    # Clean up any distortion parameters if they exist
    if wcs.sip is not None:
        wcs = wcs.celestial
    return wcs

def plot_setup(ax, wcs, title):
    """Configure common plot elements"""
    ra = ax.coords['ra']
    dec = ax.coords['dec']
    
    ra.set_axislabel('Right Ascension (ICRS)')
    dec.set_axislabel('Declination (ICRS)')
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm:ss')
    ra.set_ticks(size=10, color='white', width=1.5)
    dec.set_ticks(size=10, color='white', width=1.5)
    ax.coords.grid(True, color='white', linestyle='--', alpha=0.7)
    ax.set_title(title)
    ax.invert_xaxis()

if __name__ == '__main__':
    args = parse_args()
    target_coord = get_reference_coords(args.emodepath)

    fig = plt.figure(figsize=(20, 8))
    plt.subplots_adjust(wspace=0.4)

    # Process COADD image
    with fits.open(args.coaddpath) as hdul:
        coadd_data = hdul[0].data
        wcs_coadd = configure_wcs(hdul[0].header, target_coord)
        print(wcs_coadd)
        
        ax1 = fig.add_subplot(121, projection=wcs_coadd)
        norm1 = ImageNormalize(coadd_data, interval=ZScaleInterval())
        im1 = ax1.imshow(coadd_data, cmap='magma', norm=norm1, origin='lower')
        plot_setup(ax1, wcs_coadd, 'COADD Image')

    # Process EMODE image
    with fits.open(args.emodepath) as hdul:
        emode_data = hdul[0].data
        wcs_emode = configure_wcs(hdul[0].header, target_coord)
        
        ax2 = fig.add_subplot(122, projection=wcs_emode)
        norm2 = ImageNormalize(emode_data, interval=ZScaleInterval())
        im2 = ax2.imshow(emode_data, cmap='viridis', norm=norm2, origin='lower')
        plot_setup(ax2, wcs_emode, 'E-mode Reconstruction')
        ax2.coords['dec'].set_axislabel(' ')  # Avoid duplicate label

    # Add colorbars
    cbar1 = fig.colorbar(im1, ax=ax1, pad=0.01)
    cbar2 = fig.colorbar(im2, ax=ax2, pad=0.01)
    cbar1.set_label('Flux (Jy/beam)', fontsize=12)
    cbar2.set_label('E-mode Amplitude', fontsize=12)

    os.makedirs('./outputs', exist_ok=True)
    plt.savefig(os.path.join('./outputs', 'combined_visualization.png'), bbox_inches='tight')
    plt.close()