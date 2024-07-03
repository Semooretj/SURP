import numpy as np
import healpy as hp
import astropy.units as u
import matplotlib.pyplot as plt
import dropbox

# Needs a dropbox access token to dump maps into dropbox
dbx = dropbox.Dropbox(input('Enter your dropbox access token: '))

def compare_rotator_methods(map_index,plank_map,rot):
    """ Compare the two methods of rotating maps
    Parameters:
    map_index -------- index of the map from Plank Data, int
    plank_map -------- the Plank map data, array-like
    rot -------------- the euler rotation to apply, array-like
    
    Returns:
    map_index -------- index of the map from Plank Data, int
    ringed ----------- the original map in ringed order, array-like
    residual_alms ---- the residuals of spherical harmonics rotation, array-like
    residual_pixel --- the residuals of pixel space rotation, array-like
    """
    global dbx
    map_names = ['217 GHz', '353 GHz', '545 GHz', '857 GHz', '2998 GHz']
    
    # Convert from NESTED to RINGED
    ringed = hp.reorder(plank_map[map_index],n2r=True)
    
    # Create rotators
    rotate = hp.Rotator(rot, inv = True)  # Apparently, more accurate to use opposite to centre and invert
    rotate_back = hp.Rotator(rot)
    
    # Rotate and return the map in spherical harmonice(alms)
    shift_alms = rotate.rotate_map_alms(ringed)
    return_alms = rotate_back.rotate_map_alms(shift_alms)
    
    # Rotate and return the map in pixel space
    shift_pixel = rotate.rotate_map_pixel(ringed)
    return_pixel = rotate_back.rotate_map_pixel(shift_pixel)
    
    # Residuals
    residual_alms = return_alms - ringed
    residual_pixel = return_pixel - ringed
    
    
    # Rotated maps
    hp.mollview(np.log10(shift_alms), title=map_names[map_index]+ ' Rotated in spherical harmonics space', max = 4, min= -3 ,unit='log10(MJy/sr)')
    plt.savefig('maps/rotated_alms.pdf')
    hp.mollview(np.log10(shift_pixel), title=map_names[map_index]+' Rotated in pixel space', max = 4, min= -3 ,unit='log10(MJy/sr)')
    plt.savefig('maps/rotated_pixel.pdf')
    hp.mollview(np.log10(return_alms), title=map_names[map_index]+' Rotated back in spherical harmonics space', max = 4, min= -3 ,unit='log10(MJy/sr)')
    plt.savefig('maps/rotated_back_alms.pdf')
    hp.mollview(np.log10(return_pixel), title=map_names[map_index]+' Rotated back in pixel space', max = 4, min= -3 ,unit='log10(MJy/sr)')
    plt.savefig('maps/rotated_back_pixel.pdf')
    
    # Save to dropbox
    with open('maps/rotated_alms.pdf', 'rb') as f: im1 = f.read()
    with open('maps/rotated_pixel.pdf', 'rb') as f: im2 = f.read()
    with open('maps/rotated_back_alms.pdf', 'rb') as f: im3 = f.read()
    with open('maps/rotated_back_pixel.pdf', 'rb') as f: im4 = f.read()
    
    
    dbx.files_upload(im1, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/rotated_alms.pdf')
    dbx.files_upload(im2, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/rotated_pixel.pdf')
    dbx.files_upload(im3, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/rotated_back_alms.pdf')
    dbx.files_upload(im4, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/rotated_back_pixel.pdf')

    return map_index, ringed, residual_alms, residual_pixel

def plot_residuals_ratios(map_index,ringed_map,residual_alms,residual_pixel):
    """ Plot and save residuals and ratios
    Parameters:
    map_index -------- index of the map from Plank Data, int
    ringed_map ------- the original map in ringed order, array-like
    residual_alms ---- the residuals of spherical harmonics rotation, array-like
    residual_pixel --- the residuals of pixel space rotation, array-like
    
    Plots and saves the figures, sends them to dropbox
    """
    global dbx
    map_names = ['217 GHz', '353 GHz', '545 GHz', '857 GHz', '2998 GHz']
   
    # Plot and save alms residuals
    log_residual_alms = np.log10(np.abs(residual_alms))
    residual_alms_mask = np.isfinite(log_residual_alms)
    finite_residual_alms = log_residual_alms[residual_alms_mask]
    hp.mollview(log_residual_alms, title=map_names[map_index]+' Log10 harmonics residuals ', min = np.min(finite_residual_alms), max = np.max(finite_residual_alms) ,unit='log10(MJy/sr)')
    plt.savefig('maps/log_residuals_alms.pdf')
    with open('maps/log_residuals_alms.pdf', 'rb') as f: im3 = f.read()
    dbx.files_upload(im3, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/log_residuals_alms.pdf')
    
    # Plot and save pixel-space residuals
    log_residual_pixel = np.log10(np.abs(residual_pixel))
    residual_pixel_mask = np.isfinite(log_residual_pixel)
    finite_residual_pixel = log_residual_pixel[residual_pixel_mask]
    hp.mollview(log_residual_pixel, title=map_names[map_index]+' Log10 pixel-space residuals ', min = np.min(finite_residual_pixel), max = np.max(finite_residual_pixel) ,unit='log10(MJy/sr)')
    plt.savefig('maps/log_residuals_pixel.pdf')
    with open('maps/log_residuals_pixel.pdf', 'rb') as f: im4 = f.read()
    dbx.files_upload(im4, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/log_residuals_pixel.pdf')

    # Plot and save alms ratios
    log_ratio_alms = np.log10(np.abs(residual_alms/ringed_map))
    alms_mask = np.isfinite(log_ratio_alms)
    finite_alms = log_ratio_alms[alms_mask]
    hp.mollview(log_ratio_alms, title=map_names[map_index]+' Log10 ratio of harmonics residual to original map ', min = np.min(finite_alms), max = np.max(finite_alms))                                                                                                                          
    plt.savefig('maps/log_ratio_alms.pdf')
    with open('maps/log_ratio_alms.pdf', 'rb') as f: im5 = f.read()
    dbx.files_upload(im5, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/log_ratio_alms.pdf')
    
    # Plot and save pixel-space ratios    
    log_ratio_pixel = np.log10(np.abs(residual_pixel/ringed_map))
    pixel_mask = np.isfinite(log_ratio_pixel)
    finite_pixel = log_ratio_pixel[pixel_mask]
    hp.mollview(log_ratio_pixel, title=map_names[map_index]+' Log10 ratio of pixel-space residual to original map ', min = np.min(finite_pixel), max = np.max(finite_pixel))                                                                                                              
    plt.savefig('maps/log_ratio_pixel.pdf')
    with open('maps/log_ratio_pixel.pdf', 'rb') as f: im6 = f.read()
    dbx.files_upload(im6, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/log_ratio_pixel.pdf')
    
    return None

def ratio_histogram(map_index,ringed_map,residual_alms,residual_pixel):
    """ Plot and save histogram of ratios
    Parameters
    map_index -------- index of the map from Plank Data, int
    ringed_map ------- the original map in ringed order, array-like
    residual_alms ---- the residuals of spherical harmonics rotation, array-like
    residual_pixel --- the residuals of pixel space rotation, array-like
    
    Plots and saves the histogram of ratios, sends it to dropbox
    
    """
    global dbx
    map_names = ['217 GHz', '353 GHz', '545 GHz', '857 GHz', '2998 GHz']
    
    # Histogram of ratios
    plt.hist(np.abs(residual_alms/ringed_map), bins=100,log = True, label='Spherical harmonics', alpha = 0.9)
    plt.hist(np.abs(residual_pixel/ringed_map), bins=100, log = True, label='Pixel space', alpha = 0.7)
    plt.xlabel('Ratio of residuals to original map')
    plt.ylabel('Number of pixels')
    plt.legend()
    plt.title(map_names[map_index]+' Ratio of residuals to original map')
    plt.savefig('maps/histogram_ratios.pdf')
    plt.show()
    with open('maps/histogram_ratios.pdf', 'rb') as f: im1 = f.read()
    dbx.files_upload(im1, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/histogram_ratios.pdf')
    
   
    return None


