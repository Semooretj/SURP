import numpy as np
import healpy as hp

def rot(nside,pixel_data,rot):
    """ Rotates data based on given euler angles
    Parameters:
    nside -------- resolution of data set, int
    pixel_data --- data value for each pixel, array of float/int
    rot ---------- euler angles of rotation given in ZYX convention , list/array-like
    
    Returns:
    rotated_map -- the pixel_data set with rearranged order, array-like
    """
    # Create rotation function
    rot = hp.Rotator(rot=rot, deg=True)

    # Get array of indices to work with
    pixel_indices = []
    for index,value in enumerate(pixel_data): pixel_indices.append(index) 
    
    # Get the angular coordinates for each pixel
    theta, phi = hp.pix2ang(nside, pixel_indices)

    # Apply the rotation to these coordinates
    theta_rot, phi_rot = rot(theta, phi)

    # print(theta)
    # print(theta_rot)
    # print(phi)
    # print(phi_rot)

    # Convert the rotated coordinates back to pixel indices
    new_indices = hp.ang2pix(nside, theta_rot, phi_rot)

    # Create a new map with the pixels rearranged
    rotated_map = pixel_data[new_indices]
    return rotated_map

