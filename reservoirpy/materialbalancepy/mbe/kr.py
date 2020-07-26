import pandas as pd 
import numpy as np
from scipy.interpolate import interp1d

class kr(pd.DataFrame):
    
    def __init__(self, *args, **kwargs):
        wet_col = kwargs.pop("index", 'sw')
        wet = kwargs.pop('wet',None)
        assert wet_col in ['sw','so']
        assert isinstance(wet,(list,np.ndarray,type(None)))
        super().__init__(*args, **kwargs)
                
        # The sw must be present in the pvt class

        if wet is not None:
            wet = np.atleast_1d(wet)
            self[wet_col] = wet
            assert self[wet_col].is_monotonic_increasing or self[wet_col].is_monotonic_decreasing , "Wet phase must be increasing or decreasing"
            self.set_index(wet,inplace=True)
        elif wet_col in self.columns:
            assert self[wet_col].is_monotonic_increasing or self[wet_col].is_monotonic_decreasing , "Wet phase must be increasing or decreasing"
            self.set_index(wet_col,inplace=True)
        elif self.index.name == wet_col:
            assert self.index.is_monotonic_increasing or self[wet_col].is_monotonic_decreasing, "Wet phase must be increasing or decreasing"
    
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
                _interpolated = interp1d(self.index,self[i], bounds_error=False,fill_value='extrapolate')(p)
                int_dict[i] = _interpolated

        int_df = pd.DataFrame(int_dict, index=p)
        int_df.index.name = 'saturation'
        return int_df 
         
    @property   
    def _constructor(self):
        return kr