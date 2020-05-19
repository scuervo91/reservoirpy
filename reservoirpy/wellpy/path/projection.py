import numpy as np 
import pandas as pd  

def unit_vector(azi:(int,float)) -> np.array:
    """
    Get the unit vector2D of a given azimuth
    Input:
        azi -> (int,float) Azimuth in Degrees
    Return:
        u -> (np.ndarray) numpy array with a shape of (2,1) with the x and y components of unit vector
    """
    assert isinstance(azi,(int,float,np.ndarray))
    alpha = 90 - azi
    alpha_rad = np.deg2rad(alpha)
    x = np.cos(alpha_rad)
    y = np.sin(alpha_rad)
    p = np.array([[x,y]])
    return p

def projection_1d(x, azi, center=None):
    """
    Get the 1D projection of a series of 2D Coordinates within a given azimuth direction

    Input:
        x -> (np.ndarray) Numpy array of shape (m,2) being m the number of coordinates
        azi -> (int,float) Azimuth in Degrees
        center -(list,np.ndarray)  list or numpy array with the center 
    Return:
        u -> (np.ndarray) numpy array with a shape of (m,1)
    """
    assert isinstance(x,np.ndarray) and x.shape[1] == 2
    assert isinstance(azi,(int,float,np.ndarray))
    assert isinstance(center,(list,np.ndarray, type(None)))
    
    if isinstance(center,type(None)):
        center = x.mean(axis=0)
    else:
        center = np.atleast_1d(center)
        assert center.shape == (2,)
    #Normalize the coordinates by substracting the average coordinates

    x = x - center

    # Get the unit vector
    u = unit_vector(azi)

    # Projection over the azimuth direction
    cv = np.squeeze(np.dot(x,u.T))

    return cv, center
