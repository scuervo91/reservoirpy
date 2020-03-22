import pandas as pd
from scipy.interpolate import interp1d

class pvt(pd.DataFrame):
    _metadata = ['rs_int','rho_int','b_int','mu_int','co_int','tens_int']
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # The pressure must be present in the pvt class
        if 'pressure' not in self.columns:
            raise ValueError("Column Pressure must be present")
        if self['pressure'].is_monotonic == False:
            raise ValueError("Pressure must be increasing")  
        
        if 'rs' in self.columns: 
            self.rs_int = interp1d(self['pressure'],self['rs'])

        if 'b' in self.columns: 
            self.b_int = interp1d(self['pressure'],self['b'])

        if 'rho' in self.columns: 
            self.rho_int = interp1d(self['pressure'],self['rho'])      

        if 'mu' in self.columns: 
            self.mu_int = interp1d(self['pressure'],self['mu'])
        
        if 'co' in self.columns: 
            self.co_int = interp1d(self['pressure'],self['co'])     

        if 'ten' in self.columns: 
            self.ten_int = interp1d(self['pressure'],self['ten']) 
         
    @property   
    def _constructor(self):
        return pvt