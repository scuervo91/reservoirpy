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
        self.api = kwargs.pop("api", None)
        self.sulphur = kwargs.pop("sulphur", None)
        self.pb = kwargs.pop("pb", None)
        self.rsb = kwargs.pop("rsb", None)
        self.sg_gas = kwargs.pop("sg_gas", None)
        self.temp = kwargs.pop("temp", None)
        self.pvt = kwargs.pop('pvt',None)
        self.bg = kwargs.pop("bg", None)

#####################################################
############## Properties ###########################

    @property
    def formation(self):
        return self._formation

    @formation.setter
    def formation(self,value):
        assert isinstance(value,(str,type(None))), f'{type(value)} not accepted. Name must be str'
        self._formation = value

    @property
    def api(self):
        return self._api

    @api.setter
    def api(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be int'
        self._api = value

    @property
    def sulphur(self):
        return self._.sulphur

    @sulphur.setter
    def sulphur(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be int'
        self._sulphur = value

    @property
    def pb(self):
        return self._pb

    @pb.setter
    def pb(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be int'
        self._pb = value

    @property
    def rsb(self):
        return self._rsb

    @rsb.setter
    def rsb(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be int'
        self._rsb = value

    @property
    def sg_gas(self):
        return self._sg_gas

    @sg_gas.setter
    def sg_gas(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be int'
        self._sg_gas = value

    @property
    def temp(self):
        return self._temp

    @temp.setter
    def temp(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be int'
        self._temp = value

    @property
    def bg(self):
        return self._bg

    @bg.setter
    def bg(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be int'
        self._bg = value

    @property
    def pvt(self):
        return self._pvt

    @pvt.setter
    def pvt(self,value):
        assert isinstance(value,(oil_pvt,type(None))), f'{type(value)} not accepted. Name must be reservoirpy.pvtpy.black_oil.oil_pvt'
        self._pvt = value

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

        if (self._pb is None) & (self._rsb is None):
            raise ValueError('Either Bubble point or Gas Oil Ratio must be defined')
        elif self._pb is None:
            self._pb = pb(rs=self._rsb,temp=self._temp,sg_gas=self._sg_gas,api=self._api,
                methods=kwargs['pb'],multiple=False, correction=True)['pb'].values
            rs_cor = rs(p=p_range,pb=self._pb,temp=self._temp,api=self._api,sg_gas=self._sg_gas,
                rsb=self._rsb,multiple=False,methods=['valarde'])
        elif self._rsb is None:
            rs_cor = rs(p=p_range,pb=self._pb,temp=self._temp,api=self._api,sg_gas=self._sg_gas,
                multiple=False,methods=kwargs['rs'])

        bo_cor = bo(p=p_range,rs=rs_cor['rs'].values,pb=self._pb,temp=self._temp,api=self._api,
            sg_gas=self._sg_gas,multiple=False,methods=kwargs['bo'])
        
        co_cor = co(p=p_range,rs=rs_cor['rs'].values,pb=self._pb,temp=self._temp,api=self._api,
            sg_gas=self._sg_gas,bo=bo_cor['bo'].values,bg=self._bg,
            method_above_pb=kwargs['co']['above_pb'],method_below_pb=kwargs['co']['below_pb'])

        muo_cor = muo(p=p_range,rs=rs_cor['rs'].values,pb=self._pb,temp=self._temp,api=self._api,
            method_above_pb=kwargs['muo']['above_pb'],method_below_pb=kwargs['muo']['above_pb'],
            method_dead=kwargs['muod'])

        rho_cor = rho_oil(p=p_range,co=co_cor['co'].values,bo=bo_cor['bo'].values,rs=rs_cor['rs'].values,
            api=self._api,pb=self._pb,multiple=False,methods=kwargs['rho'])

        _pvt = pd.concat([rs_cor,bo_cor,co_cor,muo_cor,rho_cor],axis=1)
        print(_pvt)
        self._pvt=oil_pvt(_pvt.reset_index())

        return self._pvt

            









