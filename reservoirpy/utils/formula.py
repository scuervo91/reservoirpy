import numpy as np
from scipy.interpolate import interp1d

def intercept_curves(x1,y1,x2,y2, n=20):
    x1= np.atleast_1d(x1)
    x2= np.atleast_1d(x2)
    y1= np.atleast_1d(y1)
    y2= np.atleast_1d(y2)

    #Create interpolation objects for both curves
    curve_inter_1 = interp1d(x1,y1) 
    curve_inter_2 = interp1d(x2,y2)

    # Determine x range in the curves
    mins = np.array([x1.min(),x2.min()])
    maxs = np.array([x1.max(),x2.max()])

    #lower and upper range
    lower = mins.max()
    upper = maxs.min()

    x_arr = np.linspace(lower,upper,n)

    #Interpolate the entire the range for y curves
    y1_arr = curve_inter_1(x_arr)
    y2_arr = curve_inter_2(x_arr)

    # make the difference between the y curves and establish the sign
    y_sign = np.sign(y1_arr - y2_arr)

    #diferentiate the signs
    y_dif = np.diff(y_sign)
    
    #Find the values is not zero
    ix = np.where(y_dif!=0)[0]
    
    if ix.shape == (0,):
        idx = 0
        points = np.zeros((1,2))
    else:
        points = np.zeros((ix.shape[0],2))
        for i,v in enumerate(ix):
            points[i,0] = ((x_arr[v+1]-x_arr[v])/(y_sign[v+1]-y_sign[v]))*(-y_sign[v]) + x_arr[v] 
        
        points[:,1] = curve_inter_1(points[:,0])
        idx = 1
    return points, idx