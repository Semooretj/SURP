import numpy as np
import healpy as hp

def turn(nside:int,pixel_data:list) -> list:
    """hi
    """
    # calculate number of rings for given resolution
    n_rings = 4*nside -1  # determined emperically
    
    # create empty lists for each ring
    indices =[]
    for i in range(n_rings): 
        indices.append([])

    # set first pixel's colatitude as initial test
    colat0,long = hp.pix2ang(nside,0)
    test = colat0
    
    colats = []
    n = 0
    for index,value in enumerate(pixel_data):
        colat_i,long_i = hp.pix2ang(nside,index)
        if colat_i == test:
            colats.append(colat_i)
            indices[n].append(index)
        else:
            n +=1
            colats.append(colat_i)
            test = colat_i
            indices[n].append(index)
            
    #print(indices)
    # change order of each list
    new_indices = []
    for row in indices:
        np.array(row)
        new_order = np.append(row[-1],row)
        new_order = new_order[:-1]
        new_indices.append(new_order)
    #print(new_indices)
    changed_order_indices = np.concatenate(new_indices[:])
    changed_order_values = np.array([])
    for index in changed_order_indices:
        changed_order_values = np.append(changed_order_values,pixel_data[index])
        
    #print(changed_order)
    return changed_order_values
