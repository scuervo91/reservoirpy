import pandas as pd
from scipy.interpolate import interp1d

class oil_pvt(pd.DataFrame):
    _metadata = ['rs_int','rho_int','b_int','mu_int','co_int','tens_int']
    
    def __init__(self, *args, **kwargs):
        pressure = kwargs.pop("pressure", None)
        super().__init__(*args, **kwargs)
        
        # The pressure must be present in the pvt class
        if (pressure is None) & ('pressure' in self.columns):
            assert self['pressure'].is_monotonic, "Pressure must be increasing"
            pressure = self['pressure']

        assert pressure is not None, "The pressure must be set"

        self.set_index('pressure',inplace=True)
        
        if 'rs' in self.columns: 
            self.rs_int = interp1d(self.index,self['rs'])

        if 'b' in self.columns: 
            self.b_int = interp1d(self.index,self['b'])

        if 'rho' in self.columns: 
            self.rho_int = interp1d(self.index,self['rho'])      

        if 'mu' in self.columns: 
            self.mu_int = interp1d(self.index,self['mu'])
        
        if 'co' in self.columns: 
            self.co_int = interp1d(self.index,self['co'])     

        if 'ten' in self.columns: 
            self.ten_int = interp1d(self.index,self['ten']) 
         
    @property   
    def _constructor(self):
        return pvt