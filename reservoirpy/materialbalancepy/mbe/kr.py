import pandas as pd 
import numpy as np
from scipy.interpolate import interp1d

class kr(pd.DataFrame):
    
    def __init__(self, *args, **kwargs):
        sw = kwargs.pop("sw", None)
        assert isinstance(sw,(list,np.ndarray,type(None)))
        super().__init__(*args, **kwargs)
                
        # The sw must be present in the pvt class

        if sw is not None:
            sw = np.atleast_1d(sw)
            self['sw'] = sw
            assert self['sw'].is_monotonic_increasing , "sw must be increasing"
            self.set_index('sw',inplace=True)
        elif 'sw' in self.columns:
            assert self['sw'].is_monotonic_increasing , "sw must be increasing"
            self.set_index('sw',inplace=True)
        elif self.index.name == 'sw':
            assert self.index.is_monotonic_increasing, "sw must be increasing"
    
    ## Methods

    def interpolate(self,value,property=None):
        assert isinstance(value, (int,list,float,np.ndarray))
        p = np.atleast_1d(value)

        assert isinstance(property,(str,list,type(None)))

        properties = []

        if isinstance(property, str):
            properties.append(property)
        elif isinstance(property, list):
            properties.extend(property)
        else:
            properties.extend(self.columns)

        int_dict = {}

        for i in properties:
            if i in self.columns:
                _interpolated = interp1d(self.index,self[i])(p)
                int_dict[i] = _interpolated

        int_df = pd.DataFrame(int_dict, index=p)
        int_df.index.name = 'sw'
        return int_df 
         
    @property   
    def _constructor(self):
        return kr