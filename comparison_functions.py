import numpy as np
import healpy as hp
import astropy.units as u
import matplotlib.pyplot as plt
import dropbox
import os

# Needs a dropbox access token to dump maps into dropbox
dbx = dropbox.Dropbox(input('Enter your dropbox access token: '))

def compare_rotator_methods(map,rot, plank_index = None):
    """ Compare the two methods of rotating maps
    Parameters:
    map ---------------  the map data, array-like
    rot -------------- the euler rotation to apply, array-like
    plank_index -------- index of the map if from Plank Data, int from 0-4
    
    Returns:
    ringed ----------- the original map in ringed order, array-like
    residual_alms ---- the residuals of spherical harmonics rotation, array-like
    residual_pixel --- the residuals of pixel space rotation, array-like
    plank_index -------- index of the map if from Plank Data, int from 0-4
    """
    # Dropbox access token
    global dbx
  
    # Set map data and convert from NESTED to RINGED if it's a plank map
    if plank_index is None:
        ringed = map
    else:
        map = map[plank_index]
        ringed = hp.reorder(map,n2r=True)
    
       
    
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
    
      
    # To help naming if rotating plank data
    plank_map_names = ['217 GHz', '353 GHz', '545 GHz', '857 GHz', '2998 GHz']
    # Differentiate between plank data and any other maps; set naming prefix,mollview parameters, and create directories differently
    if plank_index is None:
        prefix = 'Map'
        # To ensure mollview parameters are set by default for non-plank maps
        max = None
        min = None
        unit = None
        # Create directory for non-plank maps
        if not os.path.exists('maps/Map'):
            os.makedirs('maps/Map')
        
    else: 
        # Naming for plank maps (given the index)
        prefix = plank_map_names[plank_index]
        # Set mollview parameters for plank maps
        max = 4
        min = -3
        unit = 'log10(MJy/sr)'
        # Create directory for plank maps
        if not os.path.exists('maps/'+prefix):
            os.makedirs('maps/'+prefix)

        
    # Rotated maps
    hp.mollview(np.log10(shift_alms), title= 'log10 '+ prefix + ' Rotated in spherical harmonics space', max = max, min= min ,unit=unit)
    plt.savefig('maps/'+prefix+'/rotated_alms.pdf')
    hp.mollview(np.log10(shift_pixel), title= 'log10 '+ prefix +' Rotated in pixel space', max = max, min= min ,unit=unit)
    plt.savefig('maps/'+prefix+'/rotated_pixel.pdf')
    hp.mollview(np.log10(return_alms), title= 'log10 '+ prefix +' Rotated back in spherical harmonics space', max = max, min= min ,unit=unit)
    plt.savefig('maps/'+prefix+'/rotated_back_alms.pdf')
    hp.mollview(np.log10(return_pixel), title= 'log10 '+ prefix +' Rotated back in pixel space', max = max, min= min ,unit=unit)
    plt.savefig('maps/'+prefix+'/rotated_back_pixel.pdf')
    
    # Read image files
    with open('maps/'+prefix+'/rotated_alms.pdf', 'rb') as f: im1 = f.read()
    with open('maps/'+prefix+'/rotated_pixel.pdf', 'rb') as f: im2 = f.read()
    with open('maps/'+prefix+'/rotated_back_alms.pdf', 'rb') as f: im3 = f.read()
    with open('maps/'+prefix+'/rotated_back_pixel.pdf', 'rb') as f: im4 = f.read()
    
    # Save images to dropbox
    dbx.files_upload(im1, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/'+prefix+'/rotated_alms.pdf')
    dbx.files_upload(im2, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/'+prefix+'/rotated_pixel.pdf')
    dbx.files_upload(im3, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/'+prefix+'/rotated_back_alms.pdf')
    dbx.files_upload(im4, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/'+prefix+'/rotated_back_pixel.pdf')

    return ringed, residual_alms, residual_pixel, plank_index

def plot_residuals_ratios(ringed_map,residual_alms,residual_pixel, plank_index = None):
    """ Plot and save residuals and ratios
    Parameters:
    ringed_map ------- the original map in ringed order, array-like
    residual_alms ---- the residuals of spherical harmonics rotation, array-like
    residual_pixel --- the residuals of pixel space rotation, array-like

    
    Plots and saves the figures, sends them to dropbox
    """
    # Dropbox access token
    global dbx
 
    # To help naming if rotating plank data
    plank_map_names = ['217 GHz', '353 GHz', '545 GHz', '857 GHz', '2998 GHz']
    # differentiate between plank data and any other maps; set naming prefix and create directories differently
    if plank_index is None:
        prefix = 'Map'
        unit = None # Ensure unit is set default for non-plank maps
        # Create directory for non-plank maps
        if not os.path.exists('maps/Map'):
            os.makedirs('maps/Map')
    else: 
        # Naming for plank maps (given the index)
        prefix = plank_map_names[plank_index]
        unit = 'log10(MJy/sr)'
        # Create directory for plank maps
        if not os.path.exists('maps/'+prefix):
            os.makedirs('maps/'+prefix)
   
    # Plot and save alms residuals
    log_residual_alms = np.log10(np.abs(residual_alms))
    residual_alms_mask = np.isfinite(log_residual_alms)
    finite_residual_alms = log_residual_alms[residual_alms_mask]
    hp.mollview(log_residual_alms, title= prefix+' Log10 harmonics residuals ', min = np.min(finite_residual_alms), max = np.max(finite_residual_alms) ,unit=unit)
    plt.savefig('maps/'+prefix+'/log_residuals_alms.pdf')
    with open('maps/'+prefix+'/log_residuals_alms.pdf', 'rb') as f: im3 = f.read()
    dbx.files_upload(im3, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/'+prefix+'/log_residuals_alms.pdf')
    
    # Plot and save pixel-space residuals
    log_residual_pixel = np.log10(np.abs(residual_pixel))
    residual_pixel_mask = np.isfinite(log_residual_pixel)
    finite_residual_pixel = log_residual_pixel[residual_pixel_mask]
    hp.mollview(log_residual_pixel, title= prefix+' Log10 pixel-space residuals ', min = np.min(finite_residual_pixel), max = np.max(finite_residual_pixel) ,unit=unit)
    plt.savefig('maps/'+prefix+'/log_residuals_pixel.pdf')
    with open('maps/'+prefix+'/log_residuals_pixel.pdf', 'rb') as f: im4 = f.read()
    dbx.files_upload(im4, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/'+prefix+'/log_residuals_pixel.pdf')

    # Plot and save alms ratios
    log_ratio_alms = np.log10(np.abs(residual_alms/ringed_map))
    alms_mask = np.isfinite(log_ratio_alms)
    finite_alms = log_ratio_alms[alms_mask]
    hp.mollview(log_ratio_alms, title= prefix +' Log10 ratio of harmonics residual to original map ', min = np.min(finite_alms), max = np.max(finite_alms))                                                                                                                          
    plt.savefig('maps/'+prefix+'/log_ratio_alms.pdf')
    with open('maps/'+prefix+'/log_ratio_alms.pdf', 'rb') as f: im5 = f.read()
    dbx.files_upload(im5, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/'+prefix+'/log_ratio_alms.pdf')
    
    # Plot and save pixel-space ratios    
    log_ratio_pixel = np.log10(np.abs(residual_pixel/ringed_map))
    pixel_mask = np.isfinite(log_ratio_pixel)
    finite_pixel = log_ratio_pixel[pixel_mask]
    hp.mollview(log_ratio_pixel, title= prefix +' Log10 ratio of pixel-space residual to original map ', min = np.min(finite_pixel), max = np.max(finite_pixel))                                                                                                              
    plt.savefig('maps/'+prefix+'/log_ratio_pixel.pdf')
    with open('maps/'+prefix+'/log_ratio_pixel.pdf', 'rb') as f: im6 = f.read()
    dbx.files_upload(im6, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/'+prefix+'/log_ratio_pixel.pdf')
    
    return None

def ratio_histogram(ringed_map,residual_alms,residual_pixel,plank_index = None):
    """ Plot and save histogram of ratios
    Parameters:
    ringed_map ------- the original map in ringed order, array-like
    residual_alms ---- the residuals of spherical harmonics rotation, array-like
    residual_pixel --- the residuals of pixel space rotation, array-like
    plank_index ------ index of the map if from Plank Data, int from 0-4
    
    Plots and saves the histogram of ratios, sends it to dropbox
    
    """
    # Dropbox access token
    global dbx
    
    # To help naming if rotating plank data
    plank_map_names = ['217 GHz', '353 GHz', '545 GHz', '857 GHz', '2998 GHz']
    # differentiate between plank data and any other maps; set naming prefix and create directories differently 
    if plank_index is None:
        prefix = 'Map'
        # Create directory for non-plank maps
        if not os.path.exists('maps/Map'):
            os.makedirs('maps/Map')
    else: # Naming for plank maps (given the index)
        prefix = plank_map_names[plank_index]
        # Create directory for plank maps
        if not os.path.exists('maps/'+prefix):
            os.makedirs('maps/'+prefix)
   
    
    # Histogram of ratios
    plt.hist(np.abs(residual_alms/ringed_map), bins=100, log=True, label='Spherical harmonics', alpha = 0.9)
    plt.hist(np.abs(residual_pixel/ringed_map), bins=100,log=True, label='Pixel space', alpha = 0.7)
    plt.xlabel('Ratio of residuals to original map')
    plt.ylabel('Number of pixels')
    plt.legend()
    plt.title(prefix+' Ratio of residuals to original map')
    plt.savefig('maps/'+prefix+'/histogram_ratios.pdf')
    plt.show()
    with open('maps/'+prefix+'/histogram_ratios.pdf', 'rb') as f: im1 = f.read()
    dbx.files_upload(im1, '/Apps/Overleaf/Making the Next Generation of the 3D Interstellar Medium Dust Temperature Maps/figures/'+prefix+'/histogram_ratios.pdf')
    
   
    return None


