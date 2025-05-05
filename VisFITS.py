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
    with fits.open(emode_path) as hdul:
        header = hdul[0].header
        return SkyCoord(header['CRVAL1'] * u.deg,
                       header['CRVAL2'] * u.deg,
                       frame='icrs')

def configure_wcs(original_header, target_coord):
    header = original_header.copy()
    header['CRVAL1'] = target_coord.ra.deg
    header['CRVAL2'] = target_coord.dec.deg
    wcs = WCS(header, relax=True)
    if wcs.sip is not None:
        wcs = wcs.celestial
    return wcs

def create_individual_plot(data, wcs, title, filename, cmap='magma'):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=wcs)
    
    # Plot data
    norm = ImageNormalize(data, interval=ZScaleInterval())
    im = ax.imshow(data, cmap=cmap, norm=norm, origin='lower')
    
    # Coordinate formatting
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
    
    # Colorbar
    cbar = plt.colorbar(im, pad=0.15)
    cbar.set_label('Flux (Jy/beam)' if 'COADD' in title else 'E-mode Amplitude')
    
    # Save and clean up
    os.makedirs('./outputs', exist_ok=True)
    plt.savefig(os.path.join('./outputs', filename), bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    args = parse_args()
    target_coord = get_reference_coords(args.emodepath)

    # Process and save COADD image
    with fits.open(args.coaddpath) as hdul:
        coadd_data = hdul[0].data
        wcs_coadd = configure_wcs(hdul[0].header, target_coord)
        create_individual_plot(coadd_data, wcs_coadd, 'COADD Image', 'coadd_vis.png')

    # Process and save EMODE image
    with fits.open(args.emodepath) as hdul:
        emode_data = hdul[0].data
        wcs_emode = configure_wcs(hdul[0].header, target_coord)
        create_individual_plot(emode_data, wcs_emode, 'E-mode Reconstruction', 'emode_vis.png')