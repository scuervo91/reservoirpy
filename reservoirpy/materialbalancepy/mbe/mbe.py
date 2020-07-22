import numpy as np
import pandas as pd
from ...pvtpy.black_oil import oil, water, gas 
from .kr import kr
from ...wellproductivitypy.decline import bsw_to_wor
from .aquifer import pot_aquifer
import matplotlib.pyplot as plt
import json
import os

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

def efw(m,boi,cw,swi,cf,dp):
    return (1+m)*boi*((cw*swi+cf)/(1-swi))*dp

class reservoir:
    def __init__(self,**kwargs):
        self.n = kwargs.pop('n',0)  #Original Oil in Place in bbl (barrels)
        self.g = kwargs.pop('g',0) #Original Gas in Place in  scf (Standard Cubic Feet)
        self.m = kwargs.pop('m',0) # Gas cap ratio
        self.aquifer = kwargs.pop('aquifer',None) # aquifer model
        self.oil = kwargs.pop('oil',None)
        self.water = kwargs.pop('water',None)
        self.gas = kwargs.pop('gas',None)
        self.pi = kwargs.pop('pi',0)
        self.swi = kwargs.pop('swi',0)
        self.cf = kwargs.pop('cf',0)
        self.kr_wo = kwargs.pop('kr_wo',None)
        self.k = kwargs.pop('k',0)
        self.phi = kwargs.pop('phi',0)


    # Properties
    @property
    def kr_wo(self):
        return self._kr_wo

    @kr_wo.setter 
    def kr_wo(self,value):
        assert isinstance(value,(kr,type(None)))
        self._kr_wo = value

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
        assert isinstance(value,(int,float)), "N must be numeric"
        assert value >= 0, "N must be equal o greater than 0"
        self._n = value

    @property
    def g(self):
        return self._g 

    @g.setter 
    def g(self,value):
        assert isinstance(value,(int,float)), "G must be numeric"
        assert value >= 0, "G must be equal o greater than 0"
        self._g = value

    @property
    def m(self):
        return self._m

    @m.setter 
    def m(self,value):
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

    def forecast_np(self,pressure, wp=False, winj=0, swi=None, np_max_iter=20,er_np=0.05,start_initial_conditions=True):
        """
        Make a prediction of Cummulative with a given pressure
        """
        
        # Assert pressure is One dimession
        assert isinstance(pressure,(int,float,np.ndarray))
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
        _bsw = np.zeros(pressure.shape)

        for i in range(1,len(oil_int)):
            if oil_int.index[i] >= self.oil.pb:

                #Estimate parameters from PVT table at pressure interest
                bo_p = oil_int['bo'].iloc[i]
                bo_p_minus_1 = oil_int['bo'].iloc[i-1]
                bw_p = water_int['bw'].iloc[i]
                rs_p = oil_int['rs'].iloc[i]
                rs_p_minus_1 = oil_int['rs'].iloc[i-1]
                bg_p = gas_int['bg'].iloc[i]
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

                # Numerator part of the MBE.  F = N[Eo + m*Eg + Efw] + We + Winj*Bw + Ginj*Binj
                num = (_eo + self.m*_eg + _efw) + we + water_int['winj'].iloc[i]*bw_p

                #If WOR is used 
                if _use_wor:
                    kr_int = self.kr_wo.interpolate(_sw[i-1]).iloc[0]
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
                    print(np_it)
                else:
                    _np[i] = num / bo_p

                _gp[i] = _np[i] * rs_p

                _so[i] = (1-_sw[0])*(1-_np.sum())*(bo_p/oil_int['bo'].iloc[0])
                _sw[i] = 1 - _so[i]

            _df = pd.DataFrame({'np':_np.cumsum()*self.n,'gp':_gp.cumsum()*self.n,'wp':_wp.cumsum()*self.n,'wor':_wor,'bsw':_bsw,'sw':_sw,'so':_so,'sg':_sg}, index=pressure)
        return _df



                

                


                



