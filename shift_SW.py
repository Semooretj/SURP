import numpy as np
import healpy as hp 

def shift_SW(nside,pixel_data):
    """ Shifts the entire map in north-east direction by replacing each pixel with its SW neighbour
    Parameters:
    nside -------- resolution of data set, int
    pixel_data --- data value for each pixel, array of float/int
    
    Returns:
    changed_order_data -- the pixel_data set with rearranged order
    """
    
    # Collect indices of every SW neighbour
    SW = []
    for index, value in enumerate(pixel_data):
        neighbours = hp.get_all_neighbours(nside,index)
        #print(index,neighbours)
        SW.append(neighbours[0])
    # Get values coresponding to those indices   
    changed_order_data = pixel_data[SW]
    #print(changed_order_data)
    return changed_order_data