import numpy as np
import pandas as pd 
from typing import Union

def ssrf_dp(
    q: Union[int, float, np.ndarray],
    b: Union[int, float, np.ndarray],
    mu: Union[int, float, np.ndarray],
    k: Union[int, float, np.ndarray],
    h: Union[int, float, np.ndarray],
    re: Union[int, float, np.ndarray],
    rw: Union[int, float, np.ndarray],
    s: Union[int, float, np.ndarray]    
):
    """ssrf_dp [Obtain the Pressure Drop in psi of a single phase Steady-State Radial Flow]

    Parameters
    ----------
    q : Union[int, float, np.ndarray]
        [Flow in bbl/d]
    b : Union[int, float, np.ndarray]
        [Volumetric Factor]
    mu : Union[int, float, np.ndarray]
        [Viscosity in cP]
    k : Union[int, float, np.ndarray]h
        [Permeability in md]
    k : Union[int, float, np.ndarray]h
        [height in ft]
    re : Union[int, float, np.ndarray]
        [External radius in ft]
    rw : Union[int, float, np.ndarray]
        [Wellbore Radius in ft]
    s : Union[int, float, np.ndarray]
        [Skin factor]
    """
    
    return (141.2*q*b*mu)*(1/(k*h))*(np.log(re/rw) + s)

def ssrf_q(
    dp: Union[int, float, np.ndarray],
    b: Union[int, float, np.ndarray],
    mu: Union[int, float, np.ndarray],
    k: Union[int, float, np.ndarray],
    h: Union[int, float, np.ndarray],
    re: Union[int, float, np.ndarray],
    rw: Union[int, float, np.ndarray],
    s: Union[int, float, np.ndarray]   
):
    """ssrf_dp [Obtain the flow in bbl/d of a single phase Steady-State Radial Flow]

    Parameters
    ----------
    dp : Union[int, float, np.ndarray]
        [delta pressure in psi]
    b : Union[int, float, np.ndarray]
        [Volumetric Factor]
    mu : Union[int, float, np.ndarray]
        [Viscosity in cP]
    k : Union[int, float, np.ndarray]h
        [Permeability in md]
    k : Union[int, float, np.ndarray]h
        [height in ft]
    re : Union[int, float, np.ndarray]
        [External radius in ft]
    rw : Union[int, float, np.ndarray]
        [Wellbore Radius in ft]
    s : Union[int, float, np.ndarray]
        [Skin factor]
    """
    
    return 0.00708*k*h*dp*(1/(b*mu*(np.log(re/rw)+s)))



    
