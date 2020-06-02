import pandas as pd
import numpy as np

class pot_aquifer:
    ## Pot aquifer model MBE Tarek Ahmed Reservoir Engineer Handbook
    def __init__(self,**kwargs):
        self.k = kwargs.pop('k',None)

    @property
    def k(self):
        return self._k

    @k.setter 
    def k(self, value):
        assert isinstance(value,(int, float, np.ndarray, type(None))), 'k must be a number'
        if isinstance(value, np.ndarray):
            assert value.shape == (1,), 'Shape of numpy array must be (1,)'
            value  = value.item()
        self._k = value

    def from_parameters(self, cw, cf, ra, re, h, phi, angle):
        wi = (np.pi*(np.power(ra,2)-np.power(re,2))*h*phi)/5.615
        f = angle/365
        _k = (cw + cf)*wi*f
        self._k = _k.item()

    def we(self, dp):
        assert self.k is not None
        return self.k * dp
        

