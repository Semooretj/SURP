import numpy as np
import healpy as hp
import astropy.units as u
import matplotlib.pyplot as plt
import h5py
from scipy.spatial.transform import Rotation as rot

def test_data(lower_nside, higher_nside,  location):
    """ Create test data with a superpixel highlighted at a given location
    Parameters:
    lower_nside -------- resolution of the superpixels, int
    higher_nside ------- resolution of the data pixels, int
    location: ---------- location of the superpixel that you want to mark in the higher resolution map, longitude and latitude in degrees, tuple
    
    Returns: 
    high_res_map ------- map of the higher resolution with the superpixel at the desired location marked, numpy array
    """
    # Make sure higher_nside is greater than lower_nside
    if higher_nside < lower_nside:
        raise ValueError("higher_nside must be greater than lower_nside")
    
    # Number of pixels in the map of lower_nside
    npix1 = hp.nside2npix(lower_nside)
    
    # Create an empty map for higher_nside resolution
    high_res_map = np.zeros(hp.nside2npix(higher_nside))

    # Number of pixels in each superpixel
    npix_superpixel = (higher_nside// lower_nside)**2

    # Find the index of the superpixel that is at the desired location
    superpixel_index = hp.ang2pix(int(lower_nside), location[0], location[1], lonlat=True, nest=True) 

    # Find the indices of the pixels in the superpixel
    high_res_map[superpixel_index * npix_superpixel : (superpixel_index + 1) * npix_superpixel] = 1
    hp.mollview(high_res_map, title='Superpixel at location '+str(location), nest = True, unit='MJy/sr')
    
    return high_res_map


def get_unique_positions(lower_nside,higher_nside):
    """Get the positions for the superpixel shifts in the higher resolution map
    Parameters:
    lower_nside -------- resolution of the superpixels, int
    higher_nside ------- resolution of the data pixels, int
    
    Returns:
    vec_rot ----------- the positions for the superpixel shifts, array of position arrays [long,lat]
    """
    # Get resolution in degrees
    res = hp.nside2resol(higher_nside)*u.rad.to(u.deg)
    # Number of pixels in each side of superpixel
    side_pixels = higher_nside/lower_nside
    # Set indices for pixels in the superpixel
    x = np.arange(0,side_pixels,1)
    y = np.arange(0,side_pixels,1)
    # Rotate indices to find new positions
    vec_rot = []
    for i in x:
        for j in y:
            vec_rot.append(np.array([res*(i-j)/np.sqrt(2),res*(i+j)/np.sqrt(2)]))
    
    return np.array(vec_rot)


def generate_maps(data,positions,lower_nside,higher_nside,plank = None):
    """ Create new maps as hdf5 files for each shift position
    Parameters:
    data ------------ the data array, array-like
    positions ------- the positions for the superpixel shifts, array of position arrays [long,lat]
    lower_nside ----- resolution of the superpixels, int
    higher_nside ---- resolution of the data pixels, int
    plank ----------- the plank frequency, str
    
    Returns:
    Visualises each map, saves each new map as an hdf5 file
    """
    # Convert data to ringed order to use function
    ringed_data = hp.reorder(data,n2r=True)

    for i in range(len(positions)):
        # Rotate center to each position
        rotate = hp.Rotator(positions[i],inv = True)
        rotate_pixel = rotate.rotate_map_pixel(ringed_data)    
        
        # Set file name
        if plank is None:
            filename = 'shift_'+str(int(lower_nside))+'_'+str(int(higher_nside))+'_'+str(i)+'.h5'
        else:
            filename = str(plank)+'_shift_'+str(int(lower_nside))+'_'+str(int(higher_nside))+'_'+str(i)+'.h5'
        # Create an HDF5 file
        with h5py.File('rotated_maps/'+filename, 'w') as hdf:
            # Create a dataset in the file
            hdf.create_dataset('shifted', data=rotate_pixel)
            # Add metadata
            hdf['shifted'].attrs['position'] = positions[i]
            hdf['shifted'].attrs['ordering'] = 'ringed'
            hdf['shifted'].attrs['lower_nside'] = lower_nside
            hdf['shifted'].attrs['higher_nside'] = higher_nside
            if plank is not None:
                hdf['shifted'].attrs['plank'] = plank
        # View the rotated map            
        hp.mollview(rotate_pixel, title='Rotated to position '+str(positions[i]), unit='MJy/sr')


def shift_maps(data,lower_nside,higher_nside,plank = None):
    """ Perform shifts and generate map files
    Parameters:
    data ------------ the data array, array-like
    lower_nside ----- resolution of the superpixels, int
    higher_nside ---- resolution of the data pixels, int
    plank ----------- the plank frequency, str
    
    Returns:
    Visualises each map, saves each new map as an hdf5 file
    """
    # Get shift positions
    positions = get_unique_positions(lower_nside, higher_nside)
    # Generate maps in hdf5 file
    generate_maps(data, positions, lower_nside, higher_nside,plank)
    
    
def shift_back_inv(filename):
    """ Function for shifting map back the opposite direction of original shift
    Parameters:
    filename -------- the name of the hdf5 map file, string
    
    Returns:
    Visualises the map shifted back. This method uses the negative of the original euler angles
    """
    # Read file
    with h5py.File(filename, 'r') as hdf:
        # Get the metadata
        position = hdf['shifted'].attrs['position']
        ordering = hdf['shifted'].attrs['ordering']
        # Read the data
        data = hdf['shifted'][:]
        # if ordering != 'ringed':
        #     data = hp.reorder(data,n2r=True)
        
    # Rotate center to each position
    rotate = hp.Rotator(-position,inv = True)
    rotate_pixel = rotate.rotate_map_pixel(data)
    hp.mollview(rotate_pixel, title='Rotated back with inv', unit='MJy/sr')
    return rotate_pixel
     
     
def shift_back_noinv(filename):
    """ Function for shifting map back the opposite direction of original shift
    Parameters:
    filename -------- the name of the hdf5 map file, string
    
    Returns:
    Visualises the map shifted back. This excludes the inv = True parameter from original rotation 
    """
    # Read file
    with h5py.File(filename, 'r') as hdf:
        # Get the metadata
        position = hdf['shifted'].attrs['position']
        ordering = hdf['shifted'].attrs['ordering']
        # Read the data
        data = hdf['shifted'][:]
        # if ordering != 'ringed':
        #     data = hp.reorder(data,n2r=True)
        
    # Rotate center to each position
    rotate = hp.Rotator(position)
    rotate_pixel = rotate.rotate_map_pixel(data)
    hp.mollview(rotate_pixel, title='Rotated back without inv', unit='MJy/sr')
    return rotate_pixel
    
    
    
def get_position_vectors(nside):
    """Get the position vectors of the superpixels
    Parameters:
    nside: int, the nside of the superpixels
    
    Returns:
    vecs: array, the position vectors of the superpixels
    """
    # Number of superpixels
    npix = hp.nside2npix(nside)
    # Get the current superpixel position vectors
    vecs = hp.pix2vec(nside,range(npix),nest =False)
    vecs = np.array(vecs).T
    return vecs


def unrotate_superpixels (vector, euler_angle):
    """Get angular positions of the superpixel after rotating back by the same angle
    Parameters:
    vector: array, the position vectors of the superpixels
    euler_angle: float, the angle by which the superpixels were rotated
    
    Reurns:
    angles: array, the angular positions of the superpixel after rotating back
    """   
    # Set the rotation object
    r = rot.from_euler('zy', euler_angle, degrees=True)
    # Convert to 3x3 rotation matrix
    rot_matrix = r.as_matrix()
    # Get inverse
    rot_matrix_inv = np.linalg.inv(rot_matrix)
    # Rotate the superpixels back
    vec_rotated = np.dot(rot_matrix_inv, vector.T).T
    # Convert to longitude and latitude
    angles = hp.vec2ang(vec_rotated, lonlat=True)
    angles = np.array(angles).T
    
    return angles

    
def unrotate_superpixels_file (file):
    """Get angular positions of each superpixel after rotating back by the same amount as the rotated map
    Parameters:
    file -------- the name of the hdf5 map file, string
    
    Returns:
    angles ----------- the angular positions of each superpixel after rotating back, array of position arrays [long,lat]
    """
    # Read the HDF5 file
    with h5py.File(file, 'r') as hdf:
    # Access the dataset
        for dset in hdf:
            data = hdf[dset][:]
            #hp.mollview(data, title='Test map')
            # Access metadata
            position = hdf[dset].attrs['position']
            ordering = hdf[dset].attrs['ordering']
            higher_nside = hdf[dset].attrs['higher_nside']
            lower_nside = hdf[dset].attrs['lower_nside'] 
 
    
    # Number of superpixels
    npix = hp.nside2npix(lower_nside)
    # Get the current superpixel position vectors
    vecs = hp.pix2vec(lower_nside,range(npix),nest =False)
    # List of position vectors
    vecs = np.array(vecs).T
    # Set the rotation object
    r = rot.from_euler('zy', position, degrees=True)
    # Convert to 3x3 rotation matrix
    rot_matrix = r.as_matrix()
    # Get inverse
    rot_matrix_inv = np.linalg.inv(rot_matrix)
    # Rotate the superpixels back
    vecs_rotated = np.dot(rot_matrix_inv, vecs.T).T
    # Convert to longitude and latitude
    angles = hp.vec2ang(vecs_rotated, lonlat=True)
    angles = np.array(angles).T
    
    
    return angles
