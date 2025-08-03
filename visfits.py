from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import ImageNormalize, ZScaleInterval, PercentileInterval
from astropy.wcs import WCS
from argparse import ArgumentParser
import os
from scipy.ndimage import zoom

def parse_args():
    parser = ArgumentParser(description='Create overlay visualization of coadd and E-mode FITS images')
    parser.add_argument('coaddpath', type=str, help='Full path to COADD FITS file')
    parser.add_argument('emodepath', type=str, help='Full path to E-mode FITS file')
    parser.add_argument('--output-name', type=str, default='overlay', 
                       help='Base name for output files (default: overlay)')
    parser.add_argument('--coadd-alpha', type=float, default=0.7,
                       help='Transparency for coadd layer (default: 0.7)')
    parser.add_argument('--emode-alpha', type=float, default=0.8,
                       help='Transparency for E-mode layer (default: 0.8)')
    return parser.parse_args()

def detect_hemisphere_and_orientation(coadd_header, emode_header):
    """
    Detect hemisphere and determine if coadd needs inversion.
    Southern hemisphere observations often have inverted X coordinates.
    """
    # Check declination
    coadd_dec = coadd_header.get('CRVAL2', 0)
    emode_dec = emode_header.get('CRVAL2', 0)
    
    is_southern = coadd_dec < 0 or emode_dec < 0
    
    # Additional checks for coordinate system orientation
    coadd_cd11 = coadd_header.get('CD1_1', coadd_header.get('CDELT1', 1))
    coadd_cd22 = coadd_header.get('CD2_2', coadd_header.get('CDELT2', 1))
    
    # Check if coordinate system is inverted (negative pixel scale suggests inversion)
    coord_inverted = coadd_cd11 < 0
    
    print(f"Hemisphere detection:")
    print(f"  Coadd declination: {coadd_dec:.4f}°")
    print(f"  E-mode declination: {emode_dec:.4f}°")
    print(f"  Southern hemisphere: {is_southern}")
    print(f"  CD1_1 (X pixel scale): {coadd_cd11}")
    print(f"  Coordinate system inverted: {coord_inverted}")
    
    # For southern hemisphere, we typically need to flip the coadd
    needs_coadd_flip = is_southern
    
    return is_southern, needs_coadd_flip, coord_inverted

def create_process_visualization(coadd_aligned, emode_aligned, reference_wcs, 
                               output_name, coadd_alpha, emode_alpha, is_southern, coord_inverted):
    """Create a figure showing: Coadd + Mass Map = Overlay process"""
    
    # Create figure with subplots - now 5 panels with adjusted symbol panel widths
    fig = plt.figure(figsize=(25, 6))
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 0.15, 1, 0.15, 1])  # Increased symbol panel width
    
    # Normalize data consistently across all subplots
    coadd_norm = ImageNormalize(coadd_aligned, interval=ZScaleInterval())
    emode_norm = ImageNormalize(emode_aligned, interval=PercentileInterval(99))
    
    # Panel 1: Coadd only
    ax1 = fig.add_subplot(gs[0], projection=reference_wcs)
    im1 = ax1.imshow(coadd_aligned, cmap='gray', norm=coadd_norm, origin='lower')
    ax1.set_title('Baryonic Matter\n(Coadd)', fontsize=14, fontweight='bold')
    
    # Panel 2: Plus symbol - moved to right edge
    ax2 = fig.add_subplot(gs[1])
    ax2.text(2.0, 0.5, '+', fontsize=60, ha='right', va='center', fontweight='bold')  # Right-aligned
    ax2.set_axis_off()
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    
    # Panel 3: Mass map only  
    ax3 = fig.add_subplot(gs[2], projection=reference_wcs)
    im3 = ax3.imshow(emode_aligned, cmap='plasma', norm=emode_norm, origin='lower')
    ax3.set_title('Dark Matter\n(E-mode)', fontsize=14, fontweight='bold')
    
    # Panel 4: Equals symbol - moved to right edge
    ax4 = fig.add_subplot(gs[3])
    ax4.text(1.8, 0.5, '=', fontsize=60, ha='right', va='center', fontweight='bold')  # Right-aligned
    ax4.set_axis_off()
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    
    # Panel 5: Final overlay
    ax5 = fig.add_subplot(gs[4], projection=reference_wcs)
    im5a = ax5.imshow(coadd_aligned, cmap='gray', norm=coadd_norm, 
                     origin='lower', alpha=coadd_alpha)
    im5b = ax5.imshow(emode_aligned, cmap='plasma', norm=emode_norm, 
                     origin='lower', alpha=emode_alpha)
    ax5.set_title('Combined Overlay\n(Dark + Baryonic Matter)', fontsize=14, fontweight='bold')
    
    # Apply coordinate handling to all WCS subplots
    for ax in [ax1, ax3, ax5]:
        ra = ax.coords['ra']
        dec = ax.coords['dec']
        ra.set_axislabel('RA', fontsize=10)
        dec.set_axislabel('Dec', fontsize=10)
        ra.set_major_formatter('hh:mm:ss')
        dec.set_major_formatter('dd:mm')
        
        # Handle coordinate display
        if is_southern and not coord_inverted:
            ax.invert_xaxis()
        elif coord_inverted and not is_southern:
            ax.invert_xaxis()
            
        # Grid
        ra.set_ticks(size=6, color='white', width=1)
        dec.set_ticks(size=6, color='white', width=1)
        ax.coords.grid(True, color='white', linestyle='--', alpha=0.4)
    
    # Add colorbars to relevant subplots
    cbar1 = plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    cbar1.set_label('Flux\n(Jy/beam)', fontsize=10)
    
    cbar3 = plt.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04)
    cbar3.set_label('E-mode\nAmplitude', fontsize=10)
    
    # Add legend to overlay
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='gray', alpha=coadd_alpha, label='Baryonic Matter'),
        Patch(facecolor='purple', alpha=emode_alpha, label='Dark Matter')
    ]
    ax5.legend(handles=legend_elements, loc='upper right', fontsize=10, 
               fancybox=True, framealpha=0.8)
    
    # Overall title
    fig.suptitle('Weak Lensing Analysis: Mass Distribution Overlay Process', 
                fontsize=16, fontweight='bold', y=0.95)
    
    # Adjust layout with tighter horizontal spacing
    plt.tight_layout()
    plt.subplots_adjust(top=0.75, wspace=0.02)  # Reduced wspace
    
    # Save process visualization
    process_filename = f'{output_name}_process.png'
    output_path = os.path.join('./outputs', process_filename)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    return process_filename

def align_images(coadd_data, coadd_wcs, emode_data, emode_wcs):
    """Align images to same grid, using larger image as reference"""
    coadd_shape = coadd_data.shape
    emode_shape = emode_data.shape
    
    # Use the larger image as reference
    if np.prod(coadd_shape) >= np.prod(emode_shape):
        # Resample E-mode to match coadd
        zoom_factors = [coadd_shape[i] / emode_shape[i] for i in range(2)]
        emode_aligned = zoom(emode_data, zoom_factors, order=1, mode='constant', cval=0)
        return coadd_data, emode_aligned, coadd_wcs
    else:
        # Resample coadd to match E-mode  
        zoom_factors = [emode_shape[i] / coadd_shape[i] for i in range(2)]
        coadd_aligned = zoom(coadd_data, zoom_factors, order=1, mode='constant', cval=0)
        return coadd_aligned, emode_data, emode_wcs

def main():
    args = parse_args()
    
    print(f"Loading FITS files...")
    print(f"Coadd: {args.coaddpath}")
    print(f"E-mode: {args.emodepath}")
    
    # Check if files exist
    if not os.path.exists(args.coaddpath):
        print(f"ERROR: Coadd file not found: {args.coaddpath}")
        return 1
    if not os.path.exists(args.emodepath):
        print(f"ERROR: E-mode file not found: {args.emodepath}")
        return 1
    
    # Load FITS data
    with fits.open(args.coaddpath) as hdul:
        coadd_data = hdul[0].data
        coadd_header = hdul[0].header
        coadd_wcs = WCS(coadd_header, relax=True)
        if coadd_wcs.sip is not None:
            coadd_wcs = coadd_wcs.celestial
    
    with fits.open(args.emodepath) as hdul:
        emode_data = hdul[0].data
        emode_header = hdul[0].header
        emode_wcs = WCS(emode_header, relax=True)
        if emode_wcs.sip is not None:
            emode_wcs = emode_wcs.celestial
    
    print(f"Coadd shape: {coadd_data.shape}")
    print(f"E-mode shape: {emode_data.shape}")
    
    # Detect hemisphere and orientation
    is_southern, needs_coadd_flip, coord_inverted = detect_hemisphere_and_orientation(
        coadd_header, emode_header)
    
    # Apply coadd inversion for southern hemisphere if needed
    if needs_coadd_flip:
        print("Applying X-axis flip to coadd for southern hemisphere observation...")
        coadd_data = np.fliplr(coadd_data)  # Flip left-right (X-axis)
    
    # Align images to same grid
    coadd_aligned, emode_aligned, reference_wcs = align_images(
        coadd_data, coadd_wcs, emode_data, emode_wcs)
    
    print(f"Aligned to shape: {coadd_aligned.shape}")
    
    # Print coordinate information for debugging
    print(f"\nCoordinate system information:")
    print(f"  Reference RA: {reference_wcs.wcs.crval[0]:.4f}°")
    print(f"  Reference Dec: {reference_wcs.wcs.crval[1]:.4f}°")
    try:
        print(f"  Pixel scale X: {reference_wcs.wcs.cd[0,0]*3600:.2f} arcsec/pixel")
        print(f"  Pixel scale Y: {reference_wcs.wcs.cd[1,1]*3600:.2f} arcsec/pixel")
    except:
        print(f"  Pixel scale info not available in CD matrix")
    
    # Create output directory
    os.makedirs('./outputs', exist_ok=True)
    
    # Set up the plot
    fig = plt.figure(figsize=(14, 12))
    ax = fig.add_subplot(111, projection=reference_wcs)
    
    # Normalize data for display
    coadd_norm = ImageNormalize(coadd_aligned, interval=ZScaleInterval())
    emode_norm = ImageNormalize(emode_aligned, interval=PercentileInterval(99))
    
    # Plot coadd (baryonic matter) as base layer
    im1 = ax.imshow(coadd_aligned, cmap='gray', norm=coadd_norm, 
                   origin='lower', alpha=args.coadd_alpha)
    
    # Overlay E-mode (dark matter) 
    im2 = ax.imshow(emode_aligned, cmap='plasma', norm=emode_norm, 
                   origin='lower', alpha=args.emode_alpha)
    
    # Set up coordinates
    ra = ax.coords['ra']
    dec = ax.coords['dec']
    ra.set_axislabel('Right Ascension (ICRS)', fontsize=14)
    dec.set_axislabel('Declination (ICRS)', fontsize=14)
    ra.set_major_formatter('hh:mm:ss.s')
    dec.set_major_formatter('dd:mm:ss')
    
    # Handle coordinate display based on hemisphere and inversion
    if is_southern and not coord_inverted:
        # For southern hemisphere with standard coordinates, invert display
        ax.invert_xaxis()
        print("Applied display X-axis inversion for southern hemisphere")
    elif coord_inverted and not is_southern:
        # For inverted coordinates in northern hemisphere, invert display
        ax.invert_xaxis()
        print("Applied display X-axis inversion for inverted coordinate system")
    
    # Grid and styling
    ra.set_ticks(size=10, color='white', width=2)
    dec.set_ticks(size=10, color='white', width=2)
    ax.coords.grid(True, color='white', linestyle='--', alpha=0.6)
    
    # Title and legend
    ax.set_title('Dark Matter vs Baryonic Matter Distribution', fontsize=16, pad=30)
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='gray', alpha=args.coadd_alpha, label='Baryonic Matter (Coadd)'),
        Patch(facecolor='purple', alpha=args.emode_alpha, label='Dark Matter (E-mode)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=12, 
             fancybox=True, shadow=True, framealpha=0.8)
    
    # Add colorbars
    cbar1 = plt.colorbar(im1, ax=ax, pad=0.02, fraction=0.046, aspect=30)
    cbar1.set_label('Flux (Jy/beam)', fontsize=12)
    cbar1.ax.yaxis.set_label_position('left')
    cbar1.ax.yaxis.tick_left()
    
    cbar2 = plt.colorbar(im2, ax=ax, pad=0.08, fraction=0.046, aspect=30)
    cbar2.set_label('E-mode Amplitude', fontsize=12)
    
    # Save the overlay
    overlay_filename = f'{args.output_name}_overlay.png'
    output_path = os.path.join('./outputs', overlay_filename)
    plt.savefig(output_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"\nSuccess! Overlay saved as: {output_path}")
    
    # Create process visualization
    print("Creating process visualization...")
    process_filename = create_process_visualization(
        coadd_aligned, emode_aligned, reference_wcs, args.output_name,
        args.coadd_alpha, args.emode_alpha, is_southern, coord_inverted)
    
    process_path = os.path.join('./outputs', process_filename)
    print(f"Process visualization saved as: {process_path}")
    
    print(f"\nGenerated files:")
    print(f"  📊 {overlay_filename} - Main overlay for your poster")
    print(f"  📈 {process_filename} - Process explanation figure")
    print(f"\nThe overlay shows dark matter distribution (purple/plasma) overlaid on baryonic matter (gray)")
    print(f"Perfect for demonstrating mass concentration around the galaxy cluster!")
    
    return 0

if __name__ == '__main__':
    exit(main())