import numpy as np
import healpy as hp


def turn(nside:int,pixel_data:list) -> list:
    """ Shifts eachjjjjj pixel to position of ring neighbour (left)
    Parameters:
    nside -------- resolution of data set
    pixel_data --- data value for each pixel
    
    Returns:
    changed_order_data -- the pixel_data set with rearranged order
    """
    # note1: maybe insert a test to confirm that size of pixel_data corresponds with nside resolution?
    
    # calculate number of rings for given resolution
    n_rings = 4*nside -1  # determined emperically
    
    # create empty lists for each ring
    indices =[]
    for i in range(n_rings): 
        indices.append([])

    # set first pixel's colatitude as initial test
    colat0,long0 = hp.pix2ang(nside,0)
    test = colat0
    
    # list to collect colatitudes -----> note2: can probably be removed from code
    colats = []
    # Ring counter
    n = 0
    # Checks if colatitude of pixel is same as test
    for index,value in enumerate(pixel_data):
        colat_i,long_i = hp.pix2ang(nside,index)
        if colat_i == test:   # Same ring
            colats.append(colat_i)  #--------> note2 contd, can probaly be removed from code
            indices[n].append(index)
        else: # Move to new ring, 
            n +=1
            colats.append(colat_i)  #--------> note2 contd, can probaly be removed from code
            test = colat_i   # change this colatitude to new test value for subsequent piels
            indices[n].append(index)
            
    #print(colats)       
    #print(indices)
    
    # change order of each list
    new_indices = []
    for row in indices:  # for each ring
        np.array(row)
        # Remove last element and make it first
        new_order = np.append(row[-1],row)
        new_order = new_order[:-1]
        new_indices.append(new_order)
    #print(new_indices)
    
    # Join all rings together to create new ordered dataset
    changed_order_indices = np.concatenate(new_indices[:])
    # Use index list to get values 
    changed_order_data = np.array([])
    for index in changed_order_indices:
        changed_order_data = np.append(changed_order_data,pixel_data[index])
    #print(changed_order)
    
    return changed_order_data

