import numpy as np
import healpy as hp
import astropy.units as u
import matplotlib.pyplot as plt


def get_unique_positions(lower_res,higher_res):
    """ """
    # Get resolution in degrees
    res = hp.nside2resol(higher_res)*u.rad.to(u.deg)
    # Number of pixels in each side of superpixel
    side_pixels = higher_res/lower_res
    # Set indices for pixels in the superpixel
    x = np.arange(0,side_pixels,1)
    y = np.arange(0,side_pixels,1)
    # Rotate indices to find new positions
    vec_rot = []
    for i in x:
        for j in y:
            vec_rot.append(np.array([res*(i-j)/np.sqrt(2),res*(i+j)/np.sqrt(2)]))
    
    return np.array(vec_rot)

def generate_maps(data,positions):
    """ """
    # Convert data to ringed order to use function
    ringed_data = hp.reorder(data,n2r=True)
    
    for i in range(len(positions)):
        # Rotate center to each position
        rotate = hp.Rotator(positions[i],inv = True)
        rotate_pixel = rotate.rotate_map_pixel(ringed_data)    
        hp.mollview(rotate_pixel, title='Rotated to position '+str(positions[i]), unit='MJy/sr')
        