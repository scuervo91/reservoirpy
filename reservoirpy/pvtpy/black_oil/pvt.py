import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from .correlations import pb, rs, co, muo, rho_oil, bo, z_factor, rhog, bg
import os

############################################################
############################################################
############################################################
## Oil PVT
class oil_pvt(pd.DataFrame):
    
    def __init__(self, *args, **kwargs):
        pressure = kwargs.pop("pressure", None)
        assert isinstance(pressure,(list,np.ndarray,type(None)))
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
    
    ## Methods

    def interpolate(self,value,property):
        assert isinstance(value, (int,float,np.ndarray))
        p = np.atleast_1d(value)

        assert isinstance(property,(str,list))

        properties = []

        if isinstance(property, str):
            properties.append(property)
        else:
            properties.extend(property)

        int_dict = {}

        for i in properties:
            if i in self.columns:
                _interpolated = interp1d(self.index,self[i])(p)
                int_dict[i] = _interpolated

        int_df = pd.DataFrame(int_dict, index=p)

        return int_df 
         
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
            'pb':'standing',
            'rs':'standing',
            'bo':'standing',
            'co':{
                'above_pb':'vazquez_beggs',
                'below_pb':'mccain'
            },
            'muod':'beal',
            'muo':{
                'above_pb':'beal',
                'below_pb':'beggs'
            },
            'rho':'banzer'
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

        self._pvt=oil_pvt(_pvt.reset_index())

        return self._pvt

           
############################################################
############################################################
############################################################
## Water PVT
class water_pvt(pd.DataFrame):
    
    def __init__(self, *args, **kwargs):
        pressure = kwargs.pop("pressure", None)
        assert isinstance(pressure,(list,np.ndarray,type(None)))
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
  
    def interpolate(self,value,property):
        assert isinstance(value, (int,float,np.ndarray))
        p = np.atleast_1d(value)

        assert isinstance(property,(str,list))

        properties = []

        if isinstance(property, str):
            properties.append(property)
        else:
            properties.extend(property)

        int_dict = {}

        for i in properties:
            if i in self.columns:
                _interpolated = interp1d(self.index,self[i])(p)
                int_dict[i] = _interpolated

        int_df = pd.DataFrame(int_dict, index=p)

        return int_df 
        
         
    @property   
    def _constructor(self):
        return water_pvt

class water:
    def __init__(self, **kwargs):

        self.formation = kwargs.pop('formation',None)
        self.pb = kwargs.pop("pb", None)
        self.salinity = kwargs.pop("s", None)
        self.temp = kwargs.pop("temp", None)
        self.pvt = kwargs.pop('pvt',None)

    #Properties
    @property
    def formation(self):
        return self._formation

    @formation.setter
    def formation(self,value):
        assert isinstance(value,(str,type(None))), f'{type(value)} not accepted. Name must be str'
        self._formation = value

    @property
    def pb(self):
        return self._pb

    @pb.setter
    def pb(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be numeric'
        assert value > 0, 'value must be possitive'
        self._pb = value

    @property
    def temp(self):
        return self._temp

    @temp.setter
    def temp(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be numeric'
        assert value > 0, 'value must be possitive'
        self._temp = value

    @property
    def salinity(self):
        return self._salinity

    @salinity.setter
    def salinity(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be numeric'
        assert value > 0, 'value must be possitive'
        self._salinity = value

    @property
    def pvt(self):
        return self._pvt

    @pvt.setter
    def pvt(self,value):
        assert isinstance(value,(water_pvt,type(None))), f'{type(value)} not accepted. Name must be reservoirpy.pvtpy.black_oil.water_pvt'
        self._pvt = value

            
############################################################
############################################################
############################################################
## Gas PVT
class gas_pvt(pd.DataFrame):
    
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
  
    def interpolate(self,value,property):
        assert isinstance(value, (int,float,np.ndarray))
        p = np.atleast_1d(value)

        assert isinstance(property,(str,list))

        properties = []

        if isinstance(property, str):
            properties.append(property)
        else:
            properties.extend(property)

        int_dict = {}

        for i in properties:
            if i in self.columns:
                _interpolated = interp1d(self.index,self[i])(p)
                int_dict[i] = _interpolated

        int_df = pd.DataFrame(int_dict, index=p)

        return int_df 
        
         
    @property   
    def _constructor(self):
        return gas_pvt


#upload table property list
file_dir = os.path.dirname(__file__)
components_path = os.path.join(file_dir,'components_properties.csv')

properties_df = pd.read_csv(components_path, index_col='name')

class chromatography(pd.DataFrame):

    def __init__(self, *args, **kwargs):
        mole_fraction = kwargs.pop("mole_fraction", None)
        assert isinstance(mole_fraction,(list,np.ndarray,type(None)))

        normalize = kwargs.pop('normalize',True)
        assert isinstance(normalize,bool)

        join = kwargs.pop('join',True)
        assert isinstance(join,bool)

        super().__init__(*args, **kwargs)

        # The Mole fraction must be present in the pvt class

        if mole_fraction is not None:
            mole_fraction = np.atleast_1d(mole_fraction)
            self['mole_fraction'] = mole_fraction
        else:
            assert 'mole_fraction' in self.columns 

        if normalize:
            _sum_mf = self['mole_fraction'].sum()
            self['mole_fraction'] = self['mole_fraction']/_sum_mf

        if join:
            _merged = properties_df.merge(self, how='inner', left_index=True, right_index=True)
            for i in _merged.columns:
                if i not in self.columns:
                    self[i] = _merged[i]
    @property
    def ma(self):
        return (self['mole_fraction']*self['mw']).sum()

    @property  
    def gas_sg(self):
        return (self['mole_fraction']*self['mw']).sum()/28.96

    @property  
    def ppc(self):
        return (self['mole_fraction']*self['ppc']).sum()

    @property  
    def tpc(self):
        return (self['mole_fraction']*self['tpc']).sum()

    def get_z(self,p=14.7,t=60, z_method='papay'):
        _ppc = self.ppc
        _tpc = self.tpc
        z = z_factor(p=p, t=t, ppc=_ppc, tpc=_tpc, method=z_method)
        return z

    def get_rhog(self,p=14.7,t=60, z_method='papay',rhog_method='real_gas'):
        _ma = self.ma
        if rhog_method == 'ideal_gas':
            _rhog = rhog(p=p,ma=_ma,t=t)
        elif rhog_method == 'real_gas':
            _z = self.get_z(p=p,t=t,z_method = z_method)
            _rhog = rhog(p=p,ma=_ma,z=_z.values.reshape(-1), t=t, method=rhog_method)
        return _rhog

    def get_sv(self,p=14.7,t=60, z_method='papay',rhog_method='real_gas'):
        _ma = self.ma
        if rhog_method == 'ideal_gas':
            _rhog = rhog(p=p,ma=_ma,t=t)
            _rhog['sv'] = 1/_rhog['rhog']
        elif rhog_method == 'real_gas':
            _z = self.get_z(p=p,t=t,z_method = z_method)
            _rhog = rhog(p=p,ma=_ma,z=_z.values.reshape(1), t=t, method=rhog_method)
            _rhog['sv'] = 1/_rhog['rhog']
        return _rhog[['sv']]

        
    @property   
    def _constructor(self):
        return chromatography






