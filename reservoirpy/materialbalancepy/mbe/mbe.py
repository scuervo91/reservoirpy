import numpy as np
import pandas as pd
from ...pvtpy.black_oil import oil, water, gas 
from ...krpy import water_oil_kr, gas_oil_kr
from ...wellproductivitypy.decline import bsw_to_wor
from .aquifer import pot_aquifer
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.api as sm
import json
import os
from scipy.optimize import curve_fit, root

def production_mechanisms_plot(ax=None):
    #Create the Axex
    rmax= ax or plt.gca()
    file_dir = os.path.dirname(__file__)
    rm_path = os.path.join(file_dir,'recovery_mechanisms.json')
    with open(rm_path,'r') as f:
        rm = json.loads(f.read())
    
    color=['gray', 'darkred', 'red', 'blue','green']
    linestyle=['-', '--', ':', '-','-.']

    for i,r in enumerate(rm['recovey_mechanisms']):
        x = rm['recovey_mechanisms'][r]['rf']
        y = rm['recovey_mechanisms'][r]['pr']

        rmax.plot(
            x,
            y,
            color=color[i],
            linestyle=linestyle[i],
            label = r
        )
    rmax.set_title('Recovery Mechanims Plot')
    rmax.set_ylabel('Normalized Reservior Pressure ')
    rmax.set_xlabel('Recovery Factor')
    rmax.set_ylim([0,1.2])
    rmax.set_xlim([0,0.6])
    rmax.legend()


# Linear form of material balance functions

def f(npp,bo,rp,rs,bg,wp,bw):
    return npp*(bo+(rp-rs)*bg) + wp*bw 

def eo(bo,boi,rs,rsi,bg):
    return (bo-boi) + (rsi - rs)*bg 

def eg(boi,bg,bgi):
    return boi*((bg/bgi)-1)

def efw(boi,cw,swi,cf,dp):
    return boi*((cw*swi+cf)/(1-swi))*dp

class production_history(pd.DataFrame):
    def __init__(self, *args, **kwargs):

        #Indicate pressure column name
        pressure = kwargs.pop('pressure','pressure')

        #Init the pd.Dataframe
        super().__init__(*args, **kwargs)

        #Assert cumulative production are in the columns
        #assert all([i in self.columns for i in ['np','wp','gp']])

        #Make pressure index
        if pressure in self.columns:
            self.set_index('pressure',inplace=True)
        
    @property   
    def _constructor(self):
        return production_history

class oil_reservoir:
    def __init__(self,**kwargs):
        self.n = kwargs.pop('n',None)  #Original Oil in Place in bbl (barrels)
        self.g = kwargs.pop('g',None) #Original Gas in Place in  scf (Standard Cubic Feet)
        self.m = kwargs.pop('m',None) # Gas cap ratio
        self.aquifer = kwargs.pop('aquifer',None) # aquifer model
        self.oil = kwargs.pop('oil',None)
        self.water = kwargs.pop('water',None)
        self.gas = kwargs.pop('gas',None)
        self.pi = kwargs.pop('pi',0)
        self.swi = kwargs.pop('swi',0)
        self.cf = kwargs.pop('cf',0)
        self.kr_wo = kwargs.pop('kr_wo',None)
        self.kr_go = kwargs.pop('kr_go',None)
        self.k = kwargs.pop('k',0)
        self.phi = kwargs.pop('phi',0)
        self.production_history = kwargs.pop('production_history',None)


    # Properties
    @property
    def kr_wo(self):
        return self._kr_wo

    @kr_wo.setter 
    def kr_wo(self,value):
        assert isinstance(value,(water_oil_kr,type(None)))
        self._kr_wo = value

    # Properties
    @property
    def production_history(self):
        return self._production_history

    @production_history.setter 
    def production_history(self,value):
        if value is not None:
            assert isinstance(value,production_history)
        self._production_history = value

    # Properties
    @property
    def kr_go(self):
        return self._kr_go

    @kr_go.setter 
    def kr_go(self,value):
        assert isinstance(value,(gas_oil_kr,type(None)))
        self._kr_go = value

    @property
    def swi(self):
        return self._swi

    @swi.setter 
    def swi(self,value):
        assert isinstance(value,(int,float)), "swi must be numeric"
        assert value >= 0 and value<=1, "swi must be equal o greater than 0 and less than 0"
        self._swi = value

    # Properties
    @property
    def cf(self):
        return self._cf

    @cf.setter 
    def cf(self,value):
        assert isinstance(value,(int,float)), "cf must be numeric"
        assert value >= 0, "cf must be equal o greater than 0 and less than 0"
        self._cf = value

    @property
    def pi(self):
        return self._pi

    @pi.setter 
    def pi(self,value):
        assert isinstance(value,(int,float)), "pi must be numeric"
        assert value >= 0, "pi must be equal o greater than 0"
        self._pi = value

    @property
    def n(self):
        return self._n 

    @n.setter 
    def n(self,value):
        if value is not None:
            assert isinstance(value,(int,float)), "N must be numeric"
            assert value >= 0, "N must be equal o greater than 0"
        self._n = value

    @property
    def g(self):
        return self._g 

    @g.setter 
    def g(self,value):
        if value is not None:
            assert isinstance(value,(int,float)), "G must be numeric"
            assert value >= 0, "G must be equal o greater than 0"
        self._g = value

    @property
    def m(self):
        return self._m

    @m.setter 
    def m(self,value):
        if value is not None:
            assert isinstance(value,(int,float)), "m must be numeric"
            assert value >= 0, "m must be equal o greater than 0"
        self._m = value

    @property
    def aquifer(self):
        return self._aquifer

    @aquifer.setter 
    def aquifer(self,value):
        assert isinstance(value,(pot_aquifer, type(None))), "we must be an aquifer model"
        self._aquifer = value

    @property
    def k(self):
        return self._k

    @k.setter 
    def k(self,value):
        assert isinstance(value,(int,float)), "k must be numeric"
        assert value >= 0, "k must be equal o greater than 0"
        self._k = value

    @property
    def phi(self):
        return self._phi

    @phi.setter 
    def phi(self,value):
        assert isinstance(value,(int,float)), "phi must be numeric"
        assert value >= 0, "phi must be equal o greater than 0"
        self._phi = value

    @property
    def oil(self):
        return self._oil

    @oil.setter 
    def oil(self,value):
        assert isinstance(value,(oil,type(None))), "oil must be pvtpy.black_oil.oil"
        self._oil = value

    @property
    def gas(self):
        return self._gas

    @gas.setter 
    def gas(self,value):
        assert isinstance(value,(gas,type(None))), "gas must be pvtpy.black_oil.gas"
        self._gas = value

    @property
    def water(self):
        return self._water

    @water.setter 
    def water(self,value):
        assert isinstance(value,(water,type(None))), "water must be pvtpy.black_oil.water"
        self._water = value

    def calculate_mbe_parameters(self,**kwargs):

        #Make PVT interpolations
        oil_pvt = self.oil.pvt.interpolate(self.production_history.index.values)
        gas_pvt = self.gas.pvt.interpolate(self.production_history.index.values)
        water_pvt = self.water.pvt.interpolate(self.production_history.index.values)

        pvt = pd.concat([oil_pvt,gas_pvt,water_pvt],axis=1,ignore_index=False)

        self.production_history[pvt.columns] = pvt

        if 'we' not in self.production_history.columns:
            self.production_history['we'] = 0

        #Calculate rp
        self.production_history['rp'] = self.production_history['gp']*1000/self.production_history['np']

        #Calculate MBE parameters
        self.production_history['delta_p'] = self.pi - self.production_history.index

        oil_initial_conditions = self.oil.pvt.interpolate(self.pi)
        gas_initial_conditions = self.gas.pvt.interpolate(self.pi)

        self.production_history['F'] = self.production_history.apply(
            lambda x: f(
                x['np'],
                x['bo'],
                x['rp'],
                x['rs'],
                x['bg'],
                x['wp'],
                x['bw']
            ),
            axis=1
        )
        self.production_history['Eo'] = self.production_history.apply(
            lambda x: eo(
                x['bo'],
                oil_initial_conditions['bo'].iloc[0],
                x['rs'],
                oil_initial_conditions['rs'].iloc[0],
                x['bg']
            ), axis=1
        )
        self.production_history['Eg'] = self.production_history.apply(
            lambda x: eg(
                oil_initial_conditions['bo'].iloc[0],
                x['bg'],
                gas_initial_conditions['bg'].iloc[0],
            ),axis=1
        )
        self.production_history['Efw'] = self.production_history.apply(
            lambda x: efw(
                oil_initial_conditions['bo'].iloc[0],
                x['cw'],
                self.swi,
                self.cf,
                x['delta_p']
            ),axis=1
        )

    def ho_params(self):
        oil_initial_conditions = self.oil.pvt.interpolate(self.pi)
        gas_initial_conditions = self.gas.pvt.interpolate(self.pi)
        self.production_history['F_div_Eo'] = self.production_history['F']/self.production_history['Eo']
        self.production_history['Eg_div_Eo'] = self.production_history['Eg']/self.production_history['Eo']
        self.production_history['mEg'] = self.m*self.production_history['Eg']
        self.production_history['Eo_mEg_Efw'] = self.production_history['Eo'] + self.production_history['mEg'] + self.production_history['Efw']

        #terms Havlena & Odeh 
        # #http://www.fekete.com/SAN/WebHelp/FeketeHarmony/Harmony_WebHelp/Content/HTML_Files/Reference_Material/Analysis_Method_Theory/Material_Balance_Theory.htm
        self.production_history['ho_y'] = (self.production_history['F'] - self.production_history['we']*self.production_history['bo'])/(self.production_history['Eo'] + oil_initial_conditions['bo'].iloc[0]*self.production_history['Efw'])
        self.production_history['ho_x'] = (self.production_history['Eg'] + gas_initial_conditions['bg'].iloc[0]*self.production_history['Efw'])/(self.production_history['Eo'] + oil_initial_conditions['bo'].iloc[0]*self.production_history['Efw'])


    def fit(self, fit_m=False, m=None, factor=1,**kwargs):
        assert self.production_history is not None
        oil_initial_conditions = self.oil.pvt.interpolate(self.pi)
        gas_initial_conditions = self.gas.pvt.interpolate(self.pi)
        if fit_m:
            if m is None:
                def mbe(mbe_parameters,n,m):
                    f = n*(mbe_parameters['Eo'] + m*mbe_parameters['Eg'] + (1+m)*mbe_parameters['Efw']) + mbe_parameters['we']

                    return f.values

                popt, pcov = curve_fit(mbe, self.production_history[['Eo','Eg','Efw','we']], self.production_history['F']*factor, bounds=(0, [np.inf, np.inf]),**kwargs)
                self.n = popt[0]/factor
                self.m = popt[1]
                self.g = self.m * self.n * oil_initial_conditions['bo'].iloc[0] / gas_initial_conditions['bg'].iloc[0]

                return pcov
            else:
                def mbe(mbe_parameters,n):
                    f = n*(mbe_parameters['Eo'] + m*mbe_parameters['Eg'] + (1+m)*mbe_parameters['Efw']) + mbe_parameters['we']

                    return f.values

                popt, pcov = curve_fit(mbe, self.production_history[['Eo','Eg','Efw','we']], self.production_history['F']*factor, bounds=(0, np.inf),**kwargs)
                self.n = popt[0]/factor
                self.m = m
                self.g = self.m * self.n * oil_initial_conditions['bo'].iloc[0] / gas_initial_conditions['bg'].iloc[0]

                return pcov


        else:
            def mbe(mbe_parameters,n):
                f = n*(mbe_parameters['Eo'] + mbe_parameters['Efw']) + mbe_parameters['we']

                return f.values

            popt, pcov = curve_fit(mbe, self.production_history[['Eo','Efw','we']], self.production_history['F']*factor, bounds=(0, np.inf), **kwargs)
            self.n = popt[0]/factor
            self.m = 0
            self.g = 0
            return pcov

    def drive_index(self):

        assert self.production_history is not None

        oil_initial_conditions = self.oil.pvt.interpolate(self.pi)
        gas_initial_conditions = self.gas.pvt.interpolate(self.pi)

        self.production_history['A'] = self.production_history.apply(
            lambda x: x['np']*(x['bo'] + (x['rp'] - oil_initial_conditions['rs'].iloc[0]) * x['bg']),
            axis=1
        )   
        
        #Depletion Drive
        self.production_history['DDI'] = (self.n * self.production_history['Eo']) / self.production_history['F'] 
        
        #Segregation Drive
        self.production_history['SDI'] = (self.n * self.m * self.production_history['Eo'])/self.production_history['F']

        #Water Drive
        self.production_history['WDI'] = self.production_history['we']/self.production_history['F'] 

        #Expantion Drive
        self.production_history['EDI'] = (self.n * self.production_history['Efw'])/self.production_history['F']


    def gas_cap_mn_unknowns_plot(self,ax=None,**kwargs):
        gax = ax or plt.gca()

        def_kw = {
        'color': 'black',
        'marker':'o',
        's': 100
        }    
        for (k,v) in def_kw.items():
            if k not in kwargs:
                kwargs[k]=v

        gax.scatter(self.production_history['Eg_div_Eo'],self.production_history['F_div_Eo'],**kwargs)
        gax.set_ylabel('F/Eo')
        gax.set_xlabel('Eg/Eo')
        gax.set_title('m&n Unkowns')

    def gas_cap_mn_unknowns(self):

        mod = smf.ols(formula='F_div_Eo ~ Eg_div_Eo', data=self.production_history).fit()

        return mod

    def ho_mbe(self):
        mod = smf.ols(formula='F ~ Eo_mEg_Efw', data=self.production_history).fit()
        oil_pvti=self.oil.pvt.interpolate(self.pi)
        gas_pvti=self.gas.pvt.interpolate(self.pi)
        OOIP = round(mod.params['Eo_mEg_Efw']/1e6,2)
        OGIP = round((self.m*OOIP*1e6*oil_pvti['bo'].iloc[0]/gas_pvti['bg'].iloc[0])/1e9,2)
        print(f"Original Oil In Place {OOIP} MMbbl")
        print(f"Original Gas In Place {OGIP} Bscf")

        self.n = OOIP*1e6
        self.g = OGIP*1e9

    def ho_mbe_2(self):
        mod = smf.ols(formula='ho_y ~ ho_x', data=self.production_history).fit()
        OOIP = round(mod.params['Intercept']/1e6,2)
        OGIP = round(mod.params['ho_x']/1e9,2)
        print(f"Original Oil In Place {OOIP} MMbbl")
        print(f"Original Gas In Place {OGIP} Bscf")
        self.n = OOIP.item()*1e6,
        self.g = OGIP.item()*1e9,

    def forecast_np(self,pressure, wp=False, winj=0, swi=None, np_max_iter=20,
        er_np=0.05,start_initial_conditions=True,np_guesses=[1e-4,2e-4,3e-4]):
        """
        Make a prediction of Cummulative with a given pressure
        """
        
        # Assert pressure is One dimession
        assert isinstance(pressure,(int,float,list,np.ndarray))
        pressure = np.atleast_1d(pressure)
        assert pressure.ndim==1
        
        #Add the Initial pressure if forecast start from initial Conditions
        if start_initial_conditions:
            pressure = np.append(pressure,self.pi)

        #Sort Pressure descening order
        pressure = np.sort(pressure)[::-1]
        
        # Assert all pressure are less than initial pressure
        assert np.all(pressure <= self.pi)

        assert isinstance(winj,(int,float,np.ndarray))
        if isinstance(winj,np.ndarray):
            assert winj.shape == pressure.shape
        else:
            winj = np.full(pressure.shape,winj)

        #Interest pressure condictions
        oil_int = self.oil.pvt.interpolate(pressure)
        water_int = self.water.pvt.interpolate(pressure)
        gas_int = self.gas.pvt.interpolate(pressure)
        water_int['winj'] = winj

        _use_wor = self.kr_wo is not None if wp==True else False

        _sw = self.swi if swi is None else swi
        _sw = np.zeros(pressure.shape)
        _sw[0] = self.swi if swi is None else swi
        _so = np.zeros(pressure.shape)
        _so[0]=1-_sw[0]
        _sg = np.zeros(pressure.shape)
        _np = np.zeros(pressure.shape)
        _wp = np.zeros(pressure.shape)
        _gp = np.zeros(pressure.shape)        
        _wor = np.zeros(pressure.shape)
        _gor = np.zeros(pressure.shape)
        _bsw = np.zeros(pressure.shape)

        for i in range(1,len(oil_int)):
            #Estimate parameters from PVT table at pressure interest
            bo_p = oil_int['bo'].iloc[i]
            bo_p_minus_1 = oil_int['bo'].iloc[i-1]
            bw_p = water_int['bw'].iloc[i]
            rs_p = oil_int['rs'].iloc[i]
            rs_p_minus_1 = oil_int['rs'].iloc[i-1]
            bg_p = gas_int['bg'].iloc[i]
            mug_p = gas_int['mug'].iloc[i]
            bg_p_minus_1 = gas_int['bg'].iloc[i-1]
            cw_p = water_int['cw'].iloc[i]
            dp = oil_int.index[i]-oil_int.index[i-1]
            muo_p = oil_int['muo'].iloc[i]
            muw_p = water_int['muw'].iloc[i]

            #Estimate Linear Form MBE material balance Parameters
            _eo = eo(bo_p,bo_p_minus_1,rs_p,rs_p_minus_1,bg_p)
            _eg = eg(bo_p,bg_p,bg_p_minus_1)
            _efw = efw(self.m,bo_p_minus_1,cw_p,_sw[i],self.cf,dp)

            #If aquifer model exist call the method we with parameter dp
            we = 0 if self.aquifer is None else self.aquifer.we(dp)

            if oil_int.index[i] >= self.oil.pb:

                # Numerator part of the MBE.  F = N[Eo + m*Eg + Efw] + We + Winj*Bw + Ginj*Binj
                num = (_eo + self.m*_eg + _efw) + we + water_int['winj'].iloc[i]*bw_p

                #If WOR is used 
                if _use_wor:
                    kr_int = self.kr_wo.interpolate(_sw[i-1])['krw'].iloc[0]
                    _kro = kr_int['kro']
                    _krw = kr_int['krw']
                    _bsw[i] = 1/(1+((_kro*muw_p)/(_krw*muo_p)))
                    _wor[i] = bsw_to_wor(_bsw[i])

                    # Guess np with out wor
                    np_guess = np.zeros(np_max_iter)
                    np_guess[0] = num / bo_p
                    np_it = 0
                    e_np = 0.1

                    while e_np >= er_np and np_it < np_max_iter-1:
                        wp = np.mean((_wor[i],_wor[i-1]))*(np_guess[np_it])
                        np_guess[np_it+1] = (num - wp*bw_p)/bo_p
                        
                        #Calculate error
                        e_np = np.abs(np_guess[np_it+1]-np_guess[np_it])/np_guess[np_it+1]

                        np_it+=1
                    _np[i] = np_guess[np_it]
                    _wp[i] = wp
                else:
                    _np[i] = num / bo_p
                
                #Estimate Gp
                _gp[i] = _np[i] * rs_p
                _gor[i] = _gp[i]/_np[i]

                #Estimate Saturations
                _so[i] = (1-_sw[0])*(1-_np.sum())*(bo_p/oil_int['bo'].iloc[0])
                _sw[i] = 1 - _so[i]

            else:
                lg = len(np_guesses) # Length of np_guesses
                gp_guess1 = np.zeros(lg)
                gp_guess2 = np.zeros(lg)

                for j in range(lg):

                    # Tarners Method for Pressure below Bubble Point
                    # Reservoir Engineering Handbook Tarek Ahmed 4 Ed. pg 843
                    gp_guess1[j] = (((_eo + self.m*_eg + _efw) + we + water_int['winj'].iloc[i]*bw_p - (_np.sum()+np_guesses[j])*bo_p)/bg_p) + (_np.sum() + np_guesses[j])*rs_p
                    
                    #Estimate Saturations
                    _so_guess = (1-_sw[0])*(1-_np.sum()+np_guesses[j])*(bo_p/oil_int['bo'].iloc[0])
                    
                    #Relative Permeability Ratio Krg/kro
                    kr_ratio = self.kr_go.interpolate(_so_guess)['krg_kro'].iloc[0]

                    #Instantaneus GOR
                    gor_guess = rs_p + kr_ratio*((muo_p*bo_p)/(mug_p*bg_p))

                    #Estimate Gp
                    gp_guess2[j] = _gp[i-1] + np.mean(_gor[i-1]+gor_guess)*(np_guesses[j])
                
                
                # Fit 2 lines to a linear equation to solve
                X_reg = sm.add_constant(np.array(_np.sum() + np_guesses))
                mod1 = sm.OLS(gp_guess1,X_reg).fit()
                mod2 = sm.OLS(gp_guess2,X_reg).fit()

                mod1_params = mod1.params
                mod2_params = mod2.params
                
                #Build system lilear equations
                params_stack = np.stack([mod1_params,mod2_params], axis=0)
                _a = np.stack([params_stack[:,1]*-1,np.ones(2)], axis=1)
                _b = params_stack[:,0]

                solve_np_gp = np.linalg.solve(_a,_b)

                _np[i] = solve_np_gp[0] - _np[i-1]
                _gp[i] = solve_np_gp[1] - _gp[i-1]
                _so[i] = (1-_sw[0])*(1-_np.sum())*(bo_p/oil_int['bo'].iloc[0])
                _sg[i] = 1 - _so[i] - _sw[i-1]
                _sw[i] = _sw[i-1]
                krg_kro = self.kr_go.interpolate(_so[i])['krg_kro'].iloc[0]
                _gor[i] = rs_p + krg_kro*((muo_p*bo_p)/(mug_p*bg_p))

            _df = pd.DataFrame(
                {
                    'np':_np.cumsum()*self.n,
                    'gp':_gp.cumsum()*self.n,
                    'wp':_wp.cumsum()*self.n,
                    'wor':_wor,
                    'gor':_gor,
                    'bsw':_bsw,
                    'sw':_sw,
                    'so':_so,
                    'sg':_sg}, 
                    index=pressure
            )
        return _df

class gas_reservoir:
    def __init__(self,**kwargs):
        self.g = kwargs.pop('g',0) #Original Gas in Place in  scf (Standard Cubic Feet)
        self.aquifer = kwargs.pop('aquifer',None) # aquifer model
        self.water = kwargs.pop('water',None)
        self.gas = kwargs.pop('gas',None)
        self.pi = kwargs.pop('pi',0)
        self.swi = kwargs.pop('swi',0)
        self.cf = kwargs.pop('cf',0)
        self.kr_gw = kwargs.pop('kr_gw',None)
        self.k = kwargs.pop('k',0)
        self.phi = kwargs.pop('phi',0)
                
    # Properties
    @property
    def kr_gw(self):
        return self._kr_gw

    @kr_gw.setter 
    def kr_gw(self,value):
        if value is not None:
            assert isinstance(value,gas_oil_kr)
        self._kr_gw = value

    @property
    def swi(self):
        return self._swi

    @swi.setter 
    def swi(self,value):
        assert isinstance(value,(int,float)), "swi must be numeric"
        assert value >= 0 and value<=1, "swi must be equal o greater than 0 and less than 0"
        self._swi = value

    # Properties
    @property
    def cf(self):
        return self._cf

    @cf.setter 
    def cf(self,value):
        assert isinstance(value,(int,float)), "cf must be numeric"
        assert value >= 0, "cf must be equal o greater than 0 and less than 0"
        self._cf = value

    @property
    def pi(self):
        return self._pi

    @pi.setter 
    def pi(self,value):
        assert isinstance(value,(int,float)), "pi must be numeric"
        assert value >= 0, "pi must be equal o greater than 0"
        self._pi = value

    @property
    def g(self):
        return self._g 

    @g.setter 
    def g(self,value):
        assert isinstance(value,(int,float)), "G must be numeric"
        assert value >= 0, "G must be equal o greater than 0"
        self._g = value

    @property
    def aquifer(self):
        return self._aquifer

    @aquifer.setter 
    def aquifer(self,value):
        assert isinstance(value,(pot_aquifer, type(None))), "we must be an aquifer model"
        self._aquifer = value

    @property
    def k(self):
        return self._k

    @k.setter 
    def k(self,value):
        assert isinstance(value,(int,float)), "k must be numeric"
        assert value >= 0, "k must be equal o greater than 0"
        self._k = value

    @property
    def phi(self):
        return self._phi

    @phi.setter 
    def phi(self,value):
        assert isinstance(value,(int,float)), "phi must be numeric"
        assert value >= 0, "phi must be equal o greater than 0"
        self._phi = value

    @property
    def gas(self):
        return self._gas

    @gas.setter 
    def gas(self,value):
        assert isinstance(value,(gas,type(None))), "gas must be pvtpy.black_oil.gas"
        self._gas = value

    @property
    def water(self):
        return self._water

    @water.setter 
    def water(self,value):
        assert isinstance(value,(water,type(None))), "water must be pvtpy.black_oil.water"
        self._water = value

    #Methods
    def plot(self,
        ax=None,
        pz_kw = {},
        ann_kw = {}
        ):
    
        # Default kwargs for all Lines
        def_kw = {
        'color': 'darkred',
        'linestyle':'--',
        'linewidth': 4
        }
        for (k,v) in def_kw.items():
            if k not in pz_kw:
                pz_kw[k]=v

        #Get the ax
        pzax= ax or plt.gca()

        #Get interpolated pz
        zi = self.gas.pvt.interpolate(self.pi, property = 'z')
        
        #get pz
        pz = self.pi / zi.iloc[0,0]

        #Plot
        pzax.plot([0,self.g/1e6],[pz,0],**pz_kw)
        pzax.set_ylabel('P/Z [psi]')
        pzax.set_xlabel('Gas Cum [MMscf]')
        pzax.set_title('Gas Reservoir Material Balance Plot')
        #If plot the annotations
        ann = ann_kw.pop('ann',True)

        if ann:
            ann_fontsize = ann_kw.pop('font_size',8)
            pzax.annotate(
                f"Original Gas In Place [MMscf] \n {self.g/1e6}",
                xy = (self.g/1e6,0),
                xycoords='data',
                xytext=(0, 30), 
                textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),
                bbox={'boxstyle':'round', 'fc':'0.8'},
                fontsize = ann_fontsize
            )

            pzax.annotate(
                f"Reservoir Pressure [psi] \n {self.pi}",
                xy = (0,self.pi),
                xycoords='data',
                xytext=(30, 0), 
                textcoords='offset points',
                arrowprops=dict(arrowstyle="->"),
                bbox={'boxstyle':'round', 'fc':'0.8'},
                fontsize = ann_fontsize
            )


                


                



