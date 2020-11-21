import numpy as np
import pandas as pd 
import scipy.special as sc
from typing import Union


def iarf_ei_dp(
    q: Union[int, float, np.ndarray],
    b: Union[int, float, np.ndarray],
    mu: Union[int, float, np.ndarray],
    k: Union[int, float, np.ndarray],
    h: Union[int, float, np.ndarray],
    phi:Union[int, float, np.ndarray],
    ct:Union[int, float, np.ndarray],   
    r:Union[int, float, np.ndarray],   
    t:Union[int, float, np.ndarray]
):
    """iarf_ei_dp [Obtain the Pressure drop from Infinite-Acting Radial Flow—Ei-Function Solution]

    Parameters
    ----------
    q : Union[int, float, np.ndarray]
        [flow in bbl/d]
    b : Union[int, float, np.ndarray]
        [Volumetric factor]
    mu : Union[int, float, np.ndarray]
        [viscosity in cP]
    k : Union[int, float, np.ndarray]
        [permeability in md]
    h : Union[int, float, np.ndarray]
        [height in ft]
    phi : Union[int, float, np.ndarray]
        [porosity]
    ct : Union[int, float, np.ndarray]
        [tptal compresibility in 1/psi]
    r : Union[int, float, np.ndarray]
        [Radius in ft]
    t : Union[int, float, np.ndarray]
        [time in hours]
    """
    
    return -((70.6*q*mu*b)/(k*h)) * sc.expi((-948*phi*mu*ct*np.power(r,2))/(k*t))
    
def iarf_dp(
    q: Union[int, float, np.ndarray],
    b: Union[int, float, np.ndarray],
    mu: Union[int, float, np.ndarray],
    k: Union[int, float, np.ndarray],
    h: Union[int, float, np.ndarray],
    phi:Union[int, float, np.ndarray],
    ct:Union[int, float, np.ndarray],   
    r:Union[int, float, np.ndarray],   
    t:Union[int, float, np.ndarray],
    s:Union[int, float, np.ndarray]
):
    """iarf_ei_dp [Obtain the Pressure drop from Infinite-Acting Radial Flow]

    Parameters
    ----------
    q : Union[int, float, np.ndarray]
        [flow in bbl/d]
    b : Union[int, float, np.ndarray]
        [Volumetric factor]
    mu : Union[int, float, np.ndarray]
        [viscosity in cP]
    k : Union[int, float, np.ndarray]
        [permeability in md]
    h : Union[int, float, np.ndarray]
        [height in ft]
    phi : Union[int, float, np.ndarray]
        [porosity]
    ct : Union[int, float, np.ndarray]
        [tptal compresibility in 1/psi]
    r : Union[int, float, np.ndarray]
        [Radius in ft]
    t : Union[int, float, np.ndarray]
        [time in hours]
    s : Union[int, float, np.ndarray]
        [skin]
    """
    
    return ((162.6*q*b*mu)/(k*h))*(np.log10((k*t)/(phi*mu*ct*np.power(rw,2))-3.23+0.87*s))