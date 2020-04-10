import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from .correlations import pb, rs, co, muo, rho_oil, bo

class oil_pvt(pd.DataFrame):
    _metadata = ['rs_int','rho_int','b_int','mu_int','co_int','tens_int']
    
    def __init__(self, *args, **kwargs):
        pressure = kwargs.pop("pressure", None)
        assert(pressure,(list,np.ndarray,type(None)))
        super().__init__(*args, **kwargs)
                
        # The pressure must be present in the pvt class

        if pressure is not None:
            pressure = np.atleast_1d(pressure)
            self['pressure'] = pressure
            assert self['pressure'].is_monotonic, "Pressure must be increasing"
            self.set_index('pressure',inplace=True)
        elif 'pressure' in self.columns:
            assert self['pressure'].is_monotonic, "Pressure must be increasing"
            self.set_index('pressure',inplace=True)
        elif self.index.name == 'pressure':
            assert self.index.is_monotonic, "Pressure must be increasing"
  
    def set_interpolators(self):
        if 'rs' in self.columns: 
            self.rs_int = interp1d(self.index,self['rs'])

        if 'bo' in self.columns: 
            self.b_int = interp1d(self.index,self['bo'])

        if 'rho' in self.columns: 
            self.rho_int = interp1d(self.index,self['rho'])      

        if 'muo' in self.columns: 
            self.muo_int = interp1d(self.index,self['muo'])
        
        if 'co' in self.columns: 
            self.co_int = interp1d(self.index,self['co'])     

        if 'ten' in self.columns: 
            self.ten_int = interp1d(self.index,self['ten']) 
        
         
    @property   
    def _constructor(self):
        return oil_pvt

class oil:
    def __init__(self, **kwargs):

        self.formation = kwargs.pop('formation',None)
        assert isinstance(self.formation,(str,type(None)))

        self.api = kwargs.pop("api", None)
        assert isinstance(self.api,(int,float,np.ndarray,type(None)))

        self.sulphur = kwargs.pop("sulphur", None)
        assert isinstance(self.sulphur,(int,float,np.ndarray,type(None)))

        self.pb = kwargs.pop("pb", None)
        assert isinstance(self.pb,(int,float,np.ndarray,type(None)))

        self.rsb = kwargs.pop("rsb", None)
        assert isinstance(self.rsb,(int,float,np.ndarray,type(None)))

        self.sg_gas = kwargs.pop("sg_gas", None)
        assert isinstance(self.sg_gas,(int,float,np.ndarray,type(None)))

        self.temp = kwargs.pop("temp", None)
        assert isinstance(self.temp,(int,float,np.ndarray,type(None)))

        self.pvt = kwargs.pop('pvt',None)
        assert isinstance(self.pvt,(oil_pvt,type(None)))

        self.bg = kwargs.pop("bg", None)
        assert isinstance(self.bg,(int,float,np.ndarray,type(None)))

    def pvt_from_correlations(self,start_pressure=20,end_pressure=5000,n=20,**kwargs):

        p_range=np.linspace(start_pressure,end_pressure,n)

        def_corr = {
            'pb':['standing'],
            'rs':['standing'],
            'bo':['standing'],
            'co':{
                'above_pb':['vazquez_beggs'],
                'below_pb':['mccain']
            },
            'muod':['beal'],
            'muo':{
                'above_pb':['beal'],
                'below_pb':['beggs']
            },
            'rho':['banzer']
        }

        for (k,v) in def_corr.items():
            if k not in kwargs:
                kwargs[k]=v

        if (self.pb is None) & (self.rsb is None):
            raise ValueError('Either Bubble point or Gas Oil Ratio must be defined')
        elif self.pb is None:
            self.pb = pb(rs=self.rsb,temp=self.temp,sg_gas=self.sg_gas,api=self.api,
                methods=kwargs['pb'],multiple=False, correction=True)['pb'].values
            rs_cor = rs(p=p_range,pb=self.pb,temp=self.temp,api=self.api,sg_gas=self.sg_gas,
                rsb=self.rsb,multiple=False,methods=['valarde'])
        elif self.rsb is None:
            rs_cor = rs(p=p_range,pb=self.pb,temp=self.temp,api=self.api,sg_gas=self.sg_gas,
                multiple=False,methods=kwargs['rs'])

        bo_cor = bo(p=p_range,rs=rs_cor['rs'].values,pb=self.pb,temp=self.temp,api=self.api,
            sg_gas=self.sg_gas,multiple=False,methods=kwargs['bo'])
        
        co_cor = co(p=p_range,rs=rs_cor['rs'].values,pb=self.pb,temp=self.temp,api=self.api,
            sg_gas=self.sg_gas,bo=bo_cor['bo'].values,bg=self.bg,
            method_above_pb=kwargs['co']['above_pb'],method_below_pb=kwargs['co']['below_pb'])

        muo_cor = muo(p=p_range,rs=rs_cor['rs'].values,pb=self.pb,temp=self.temp,api=self.api,
            method_above_pb=kwargs['muo']['above_pb'],method_below_pb=kwargs['muo']['above_pb'],
            method_dead=kwargs['muod'])

        rho_cor = rho_oil(p=p_range,co=co_cor['co'].values,bo=bo_cor['bo'].values,rs=rs_cor['rs'].values,
            api=self.api,pb=self.pb,multiple=False,methods=kwargs['rho'])

        _pvt = pd.concat([rs_cor,bo_cor,co_cor,muo_cor,rho_cor],axis=1)
        print(_pvt)
        self.pvt=oil_pvt(_pvt.reset_index())

        return self.pvt

            









