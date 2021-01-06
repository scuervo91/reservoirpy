import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from .correlations import *
import os

############################################################
############################################################
############################################################
## Oil PVT
class Pvt(pd.DataFrame):
    
    def __init__(self, *args, **kwargs):
        pressure = kwargs.pop("pressure", None)
        assert isinstance(pressure,(list,np.ndarray,type(None)))
        super().__init__(*args, **kwargs)
                
        # The pressure must be present in the pvt class

        if pressure is not None:
            pressure = np.atleast_1d(pressure)
            self['pressure'] = pressure
            assert self['pressure'].is_monotonic_increasing or self['pressure'].is_monotonic_decreasing , "Pressure must be increasing"
            self.set_index('pressure',inplace=True)
        elif 'pressure' in self.columns:
            assert self['pressure'].is_monotonic_increasing or self['pressure'].is_monotonic_decreasing, "Pressure must be increasing"
            self.set_index('pressure',inplace=True)
        elif self.index.name == 'pressure':
            assert self.index.is_monotonic_increasing or self.index.is_monotonic_decreasing, "Pressure must be increasing"
    
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
                _interpolated = interp1d(self.index,self[i],bounds_error=False,fill_value='extrapolate')(p)
                int_dict[i] = _interpolated

        int_df = pd.DataFrame(int_dict, index=p)
        int_df.index.name = 'pressure'
        return int_df 
         
    @property   
    def _constructor(self):
        return Pvt

#Default Correlations
oil_def_corr = {
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

class Oil:

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
        self.correlations = kwargs.pop('correlations',oil_def_corr.copy())

    #####################################################
    ############## Properties ###########################

    @property
    def formation(self):
        return self._formation

    @formation.setter
    def formation(self,value):
        if value is not None:
            assert isinstance(value,str), f'{type(value)} not accepted. Name must be str'
        self._formation = value

    @property
    def api(self):
        return self._api

    @api.setter
    def api(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be int'
        self._api = value

    @property
    def sulphur(self):
        return self._sulphur

    @sulphur.setter
    def sulphur(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be int'
        self._sulphur = value

    @property
    def pb(self):
        return self._pb

    @pb.setter
    def pb(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be int'
        self._pb = value

    @property
    def rsb(self):
        return self._rsb

    @rsb.setter
    def rsb(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be int'
        self._rsb = value

    @property
    def sg_gas(self):
        return self._sg_gas

    @sg_gas.setter
    def sg_gas(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be int'
        self._sg_gas = value

    @property
    def temp(self):
        return self._temp

    @temp.setter
    def temp(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be int'
        self._temp = value

    @property
    def bg(self):
        return self._bg

    @bg.setter
    def bg(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be int'
        self._bg = value

    @property
    def pvt(self):
        return self._pvt

    @pvt.setter
    def pvt(self,value):
        if value is not None:
            assert isinstance(value,Pvt), f'{type(value)} not accepted. Name must be reservoirpy.pvtpy.black_oil.pvt'
        self._pvt = value

    @property
    def correlations(self):
        return self._correlations

    @correlations.setter
    def correlations(self,value):
        assert isinstance(value,dict), f'{type(value)} not accepted. Name must be Dictionary'
        self._correlations = value

    def pvt_from_correlations(self,start_pressure=20,end_pressure=5000,n=20,**kwargs):

        p_range=np.linspace(start_pressure,end_pressure,n)

        for (k,v) in self.correlations.items():
            if k not in kwargs:
                kwargs[k]=v

        if (self._pb is None) & (self._rsb is None):
            raise ValueError('Either Bubble point or Gas Oil Ratio must be defined')
        elif self._pb is None:
            self._pb = pb(rs=self._rsb,temp=self._temp,sg_gas=self._sg_gas,api=self._api,
                methods=kwargs['pb'], correction=True)['pb'].values
            rs_cor = rs(p=p_range,pb=self._pb,temp=self._temp,api=self._api,sg_gas=self._sg_gas,
                rsb=self._rsb,method='valarde')
        else:
            rs_cor = rs(p=p_range,pb=self._pb,temp=self._temp,api=self._api,sg_gas=self._sg_gas,
                methods=kwargs['rs'])

        bo_cor = bo(p=p_range,rs=rs_cor['rs'].values,pb=self._pb,temp=self._temp,api=self._api,
            sg_gas=self._sg_gas,methods=kwargs['bo'])
        
        co_cor = co(p=p_range,rs=rs_cor['rs'].values,pb=self._pb,temp=self._temp,api=self._api,
            sg_gas=self._sg_gas,bo=bo_cor['bo'].values,bg=self._bg,
            method_above_pb=kwargs['co']['above_pb'],method_below_pb=kwargs['co']['below_pb'])

        muo_cor = muo(p=p_range,rs=rs_cor['rs'].values,pb=self._pb,temp=self._temp,api=self._api,
            method_above_pb=kwargs['muo']['above_pb'],method_below_pb=kwargs['muo']['below_pb'],
            method_dead=kwargs['muod'])

        rho_cor = rho_oil(p=p_range,co=co_cor['co'].values,bo=bo_cor['bo'].values,rs=rs_cor['rs'].values,
            api=self._api,pb=self._pb,methods=kwargs['rho'])

        _pvt = pd.concat([rs_cor,bo_cor,co_cor,muo_cor,rho_cor],axis=1)

        self._pvt=Pvt(_pvt.reset_index())

        return self._pvt

    def to_ecl(
        self,
        pressure=None,
        n_sat = 10, 
        n_unsat=5, 
        min_pressure=None, 
        max_pressure=None,
        properties = ['rs','bo','muo'],
        float_format='{:.3f}'.format
    ):
        
        assert self.pvt is not None, 'PVT not defined'
        assert self.pb is not None, 'Bublle pressure not defined'
        assert self.rsb is not None, 'Rsb not defined'
        
        string = ""
        string += "-- OIL PVT TABLE FOR LIVE OIL\n"
        string += 'PVTO\n'
        string += "-- rs      pres  bo      visc\n"
        string += "-- Mscf/rb psi   RB/STB  cP  \n"
        string += "-- ------- ----  ----    ---- \n"
        
        if pressure is None:
        
            if min_pressure is None:
                min_pressure = self.pvt.index.min()
            
            if max_pressure is None:
                max_pressure = self.pvt.index.max()
            
            if min_pressure >= self.pb:
                pressure = np.linspace(min_pressure,max_pressure,n_sat)
                flag = 'unsat'
            elif max_pressure <= self.pb:
                pressure = np.linspace(min_pressure,max_pressure,n_sat)
                flag = 'sat'
            else:
                sat_pressure = np.linspace(min_pressure,self.pb,n_sat)
                unsat_pressure = np.linspace(self.pb,max_pressure, n_unsat+1)
                pressure = np.concatenate((sat_pressure,unsat_pressure[1:]))
                flag = 'mixed'
                
        pvt_df = self.pvt.interpolate(pressure,property=properties).reset_index()
        
        #convert rs units from scf/bbl to Mscf/bbl
        pvt_df['rs'] = pvt_df['rs'] * 1e-3
        
        # Write the string
        
        if flag == 'unsat':
            string += pvt_df[['rs','pressure','bo','muo']].to_string(header=False, index=False,float_format=float_format) +'\n /\n'
            
        elif flag == 'sat':
            
            for i,r in pvt_df.iterrows():
                string += pvt_df.loc[[i],['rs','pressure','bo','muo']].to_string(index=False, header=False,float_format=float_format) + '/\n'
                
            string += '/\n'
        else:
            
            #Print Saturated data
            for i,r in pvt_df[pvt_df['pressure']<self.pb].iterrows():
                string += pvt_df.loc[[i],['rs','pressure','bo','muo']].to_string(index=False, header=False,float_format=float_format) + '/\n'  
                    
            #Print data at bubble point
            string += pvt_df.loc[pvt_df['pressure']==self.pb,['rs','pressure','bo','muo']].to_string(index=False, header=False,float_format=float_format) + '\n'  
            
            string += '-- Unsaturated Data\n'
            string += pvt_df.loc[pvt_df['pressure']>self.pb,['pressure','bo','muo']].to_string(index=False, header=False,float_format=float_format)
            string += '/\n/'
        
        return string
           
############################################################
############################################################
############################################################
water_def_corr = {
    'rsw':'culberson',
    'cw': 'standing',
    'bw': 'mccain',
    'rhow':'banzer',
    'muw' : 'van_wingen'
}

class Water:
    def __init__(self, **kwargs):

        self.formation = kwargs.pop('formation',None)
        self.pb = kwargs.pop("pb", None)
        self.salinity = kwargs.pop("salinity", 0)
        self.temp = kwargs.pop("temp", None)
        self.pvt = kwargs.pop('pvt',None)
        self.correlations = kwargs.pop('correlations',water_def_corr.copy())

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
        assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be numeric'
        assert value >= 0, 'value must be possitive'
        self._salinity = value

    @property
    def pvt(self):
        return self._pvt

    @pvt.setter
    def pvt(self,value):
        assert isinstance(value,(Pvt,type(None))), f'{type(value)} not accepted. Name must be reservoirpy.pvtpy.black_oil.pvt'
        self._pvt = value

    @property
    def correlations(self):
        return self._correlations

    @correlations.setter
    def correlations(self,value):
        assert isinstance(value,dict), f'{type(value)} not accepted. Name must be Dictionary'
        self._correlations = value

    def pvt_from_correlations(self,start_pressure=20,end_pressure=5000,n=20,**kwargs):

        p_range=np.linspace(start_pressure,end_pressure,n)

        for (k,v) in self.correlations.items():
            if k not in kwargs:
                kwargs[k]=v

        rsw_cor = rsw(p=p_range, t=self.temp, s=self.salinity, method=self.correlations['rsw'])
        cw_cor = cw(p=p_range, t=self.temp, rsw=rsw_cor['rsw'].values, s=self.salinity, method=self.correlations['cw'])
        bw_cor = bw(p=p_range, t=self.temp, pb=self.pb, cw=cw_cor['cw'].values, s=self.salinity, method=self.correlations['bw'])
        rhow_cor = rhow(p=p_range,s=self.salinity, bw=bw_cor['bw'].values, method = self.correlations['rhow'])
        muw_cor = muw(p=p_range, t=self.temp, s = self.salinity,  method = self.correlations['muw'])

        _pvt = pd.concat([rsw_cor,cw_cor,bw_cor,muw_cor,rhow_cor],axis=1)

        self._pvt=Pvt(_pvt.reset_index())
############################################################
############################################################
############################################################

#upload table property list
file_dir = os.path.dirname(__file__)
components_path = os.path.join(file_dir,'components_properties.csv')

properties_df = pd.read_csv(components_path, index_col='name')

class Chromatography(pd.DataFrame):

    def __init__(self, *args, **kwargs):
        mole_fraction = kwargs.pop("mole_fraction", None)
        compound = kwargs.pop("compound", None)
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

        if compound is not None:
            assert isinstance(compound,list)
            self['compound'] = compound
            self.set_index('compound',inplace=True)
        elif 'compound' in self.columns:
            self.set_index('compound',inplace=True)

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
     
    def get_pseudo_critical_properties(self, correct=True,correct_method='wichert-aziz'):
        _ppc = (self['mole_fraction']*self['ppc']).sum()
        _tpc = (self['mole_fraction']*self['tpc']).sum()
        
        if correct:
            _co2 = self.loc['carbon-dioxide','mole_fraction'] if 'carbon-dioxide' in self.index else 0
            _n2 = self.loc['nitrogen','mole_fraction'] if 'nitrogen' in self.index else 0
            _h2s = self.loc['hydrogen-sulfide','mole_fraction'] if 'hydrogen-sulfide' in self.index else 0
            cp_dict =  critical_properties_correction(ppc=_ppc, tpc=_tpc, co2=_co2, n2=_n2, h2s=_h2s, method=correct_method)
        else:
            cp_dict = {'ppc':_ppc,'tpc':_tpc}

        return cp_dict

    def get_z(self,p=14.7,t=60, z_method='papay', cp_correction_method='wichert-aziz'):
        cp = self.get_pseudo_critical_properties(correct_method=cp_correction_method)
        z = z_factor(p=p, t=t, ppc=cp['ppc'], tpc=cp['tpc'], method=z_method)
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
        return Chromatography

gas_def_corr = {
    'cp_correction': 'wichert-aziz',
    'z':'papay',
    'rhog':'real_gas',
    'bg':{'unit':'bbl/scf'},
    'mug':'lee_gonzalez',
    'cg':'ideal_gas'
}
class Gas:
    def __init__(self, **kwargs):

        self.formation = kwargs.pop('formation',None)
        self.temp = kwargs.pop("temp", None)
        self.gas_type = kwargs.pop("gas_type",'natural_gas')
        self.pvt = kwargs.pop('pvt',None)
        self.chromatography = kwargs.pop('chromatography',None)
        self.sg = kwargs.pop('sg',None)
        self.ma = kwargs.pop('ma',None)
        self.co2 = kwargs.pop('co2',0)
        self.h2s = kwargs.pop('h2s',0)
        self.n2 = kwargs.pop('n2',0)
        self.correlations = kwargs.pop('correlations',gas_def_corr.copy())


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
    def gas_type(self):
        return self._gas_type

    @gas_type.setter
    def gas_type(self,value):
        assert isinstance(value,(str,type(None))), f'{type(value)} not accepted. Name must be str'
        assert value in ['natural_gas','condesate_gas']
        self._gas_type = value

    @property
    def temp(self):
        return self._temp

    @temp.setter
    def temp(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be numeric'
        assert value > 0, 'value must be possitive'
        self._temp = value

    @property
    def pvt(self):
        return self._pvt
    
    @pvt.setter
    def pvt(self,value):
        assert isinstance(value,(Pvt,type(None))), 'PVT must be a instance of reservoirpy.pvtpy.black_oil.pvt object'
        self._pvt = value 

    @property
    def chromatography(self):
        return self._chromatography
    
    @chromatography.setter
    def chromatography(self,value):
        assert isinstance(value,(Chromatography,type(None))), 'chromatography must be a instance of reservoirpy.pvtpy.black_oil.chromatography object'
        self._chromatography = value 

    @property
    def sg(self):
        if self.chromatography is not None:
            _sg = self.chromatography.gas_sg
        elif self._ma is not None:
            _sg = self._ma/28.96
        else:
            _sg = self._sg
        return _sg
    
    @sg.setter
    def sg(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), 'sg must be numeric'
        self._sg = value 


    @property
    def ma(self):
        if self.chromatography is not None:
            _ma = self.chromatography.ma
        elif self._sg is not None:
            _ma = self._sg*28.96
        else:
            _ma = self._ma
        return _ma
    
    @ma.setter
    def ma(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), 'ma must be numeric'
        self._ma = value 

    @property
    def co2(self):
        return self._co2 

    @co2.setter
    def co2(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), 'co2 must be numeric'
        assert value >=0 and value<=1, 'mole fraction between 0 and 1'
        self._co2 = value

    @property
    def h2s(self):
        return self._h2s

    @h2s.setter
    def h2s(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), 'h2s must be numeric'
        assert value >=0 and value<=1, 'mole fraction between 0 and 1'
        self._h2s = value

    @property
    def n2(self):
        return self._n2

    @n2.setter
    def n2(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), 'n2 must be numeric'
        assert value >=0 and value<=1, 'mole fraction between 0 and 1'
        self._n2 = value

    @property
    def correlations(self):
        return self._correlations

    @correlations.setter
    def correlations(self,value):
        assert isinstance(value,dict), f'{type(value)} not accepted. Name must be Dictionary'
        self._correlations = value

    def pseudo_critical_properties(self,correct=True,correct_method='wichert-aziz'):
        if self.chromatography is not None:
            _cp = self.chromatography.get_pseudo_critical_properties(correct=correct,correct_method=correct_method)
        elif self.sg is not None:
            _cp = critical_properties(sg=self.sg, gas_type=self.gas_type,method='standing')
            if correct:
                _cp = critical_properties_correction(ppc=_cp['ppc'], tpc=_cp['tpc'], co2=self.co2, n2=self.n2, h2s=self.h2s, method=correct_method)
        else:
            raise ValueError('Neither chromatography nor sg gas been set')
        return _cp
    

    def pvt_from_correlations(self,start_pressure=20,end_pressure=5000,n=20,**kwargs):

        p_range=np.linspace(start_pressure,end_pressure,n)

        for (k,v) in self.correlations.items():
            if k not in kwargs:
                kwargs[k]=v

        # Define Pseudo critical properties
        pcp = self.pseudo_critical_properties(self,correct_method=kwargs['cp_correction'])  #Pseudo critical properties

        # Compressibility factor z
        z_cor = z_factor(p=p_range, t=self.temp, ppc=pcp['ppc'], tpc=pcp['tpc'], method=kwargs['z'])

        # Density 
        rhog_cor = rhog(p=p_range, ma=self.ma, z=z_cor['z'].values, r=10.73, t=self.temp, method=kwargs['rhog'])
    
        #Gas volumetric factor
        bg_cor = bg(p=p_range, t=self.temp, z=z_cor['z'].values, unit=kwargs['bg']['unit'])

        #Gas viscosity
        mug_cor = mug(p=p_range,t=self.temp, rhog=rhog_cor['rhog'].values, ma=self.ma, method=kwargs['mug'])

        #Gas compressibility 
        cg_cor = cg(p=p_range, z=z_cor['z'].values, method=kwargs['cg'])

        _pvt = pd.concat([z_cor,rhog_cor,bg_cor,mug_cor,cg_cor],axis=1)

        self._pvt=Pvt(_pvt.reset_index())


    def to_ecl(
        self,
        pressure=None,
        n = 10, 
        min_pressure=None, 
        max_pressure=None,
        properties = ['bg','mug'],
        float_format='{:.3f}'.format
    ):
        
        assert self.pvt is not None, 'PVT not defined'
        
        string = ""
        string += "-- GAS PVT TABLE FOR LIVE OIL\n"
        string += 'PVDG\n'
        string += "-- pres   bg       vic  \n"
        string += "-- psi    Rb/Mscf  cP  \n"
        string += "-- ----   ----     ---- \n"
        
        if pressure is None:
        
            if min_pressure is None:
                min_pressure = self.pvt.index.min()
            
            if max_pressure is None:
                max_pressure = self.pvt.index.max()
            
            pressure = np.linspace(min_pressure,max_pressure,n)
                
        pvt_df = self.pvt.interpolate(pressure,property=properties).reset_index()
        
        #convert bo units from rb/scf to rb/Mscf
        pvt_df['bg'] = pvt_df['bg'] * 1e3
        
        # Write the string
        string += pvt_df[['pressure','bg','mug']].to_string(header=False, index=False,float_format=float_format) +'/\n'
                   
        return string  