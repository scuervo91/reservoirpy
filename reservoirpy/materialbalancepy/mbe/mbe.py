import numpy as np
import pandas as pd
from ...pvtpy.black_oil import oil, water, gas 
from .kr import kr
from ...wellproductivitypy.decline import bsw_to_wor

# Linear form of material balance functions

def f(npp,bo,rp,rs,bg,wp,bw):
    return npp*(bo+(rp-rs)*bg) + wp*bw 

def eo(bo,boi,rsi,rs,bg):
    return (bo-boi) + (rsi - rs)*bg 

def eg(boi,bg,bgi):
    return boi*((bg/bgi)-1)

def efw(m,boi,cw,swi,cf,dp):
    return (1+m)*boi*((cw*swi+cf)/(1-swi))*dp

class reservoir:
    def __init__(self,**kwargs):
        self.n = kwargs.pop('n',0)  #Original Oil in Place in MMbbl (Millions of barrels)
        self.g = kwargs.pop('g',0) #Original Gas in Place in  Bscf (Billions of Standard Cubic Feet)
        self.m = kwargs.pop('m',0) # Gas cap ratio
        self.we = kwargs.pop('we',0)
        self.wing = kwargs.pop('wing',0)
        self.ging = kwargs.pop('ging',0)
        self.oil = kwargs.pop('oil',None)
        self.water = kwargs.pop('water',None)
        self.gas = kwargs.pop('gas',None)
        self.pi = kwargs.pop('pi',0)
        self.swi = kwargs.pop('swi',0)
        self.cf = kwargs.pop('cf',0)
        self.kr_wo = kwargs.pop('kr_wo',None)


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
    def we(self):
        return self._we

    @we.setter 
    def we(self,value):
        assert isinstance(value,(int,float)), "we must be numeric"
        assert value >= 0, "we must be equal o greater than 0"
        self._we = value

    @property
    def wing(self):
        return self._wing

    @wing.setter 
    def wing(self,value):
        assert isinstance(value,(int,float)), "wing must be numeric"
        assert value >= 0, "wing must be equal o greater than 0"
        self._wing = value

    @property
    def ging(self):
        return self._ging

    @ging.setter 
    def ging(self,value):
        assert isinstance(value,(int,float)), "ging must be numeric"
        assert value >= 0, "ging must be equal o greater than 0"
        self._ging = value

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

    def forecast_np(self,pressure, wp=False):
        """
        Make a prediction of Cummulative with a given pressure
        """
        assert isinstance(pressure,(int,float,np.ndarray))
        pressure = np.atleast_1d(pressure)

        # Initial conditions pvt in pd.Series
        ic = self.oil.pvt.interpolate(self.pi).iloc[0]

        #Interest pressure condictions
        oil_int = self.oil.pvt.interpolate(pressure)
        water_int = self.water.pvt.interpolate(pressure)

        _use_wor = self.kr_wo is not None if wp==True else False

        _sw = self.swi
        forecast = pd.DataFrame()
        for i, r in oil_int.iterrows():
            if i >= self.oil.pb:
                _eo = eo(r['bo'],ic['bo'],ic['rs'],r['rs'],r['bg'])
                _eg = eg(ic['bo'],r['bg'],ic['bg'])
                dp =self.pi - i
                _efw = efw(self.m,ic['bo'],water_int.loc[i,'cw'],self.swi,self.cf,dp)

                # Numerator part of the MBE.  F = N[Eo + m*Eg + Efw] + We + Winj*Bw + Ginj*Binj
                num = self.n*(_eo + self.m*_eg + _efw) + self.we + self.wing*water_int.loc[i,'bw'] + self.ging*r['bg']
                
                #If WOR is used 
                if _use_wor:
                    kr_int = self.kr_wo.interpolate(_sw).iloc[0]
                    _kro = kr_int['kro']
                    _krw = kr_int['krw']

                    _muw = water_int.loc[i,'muw']
                    _muo = r['muo']
                    fw = 1/(1+((_kro*_muw)/(_krw*_muo)))
                    wor = bsw_to_wor(fw)
                else:
                    wor = 0
                
                
                oil_cum = num / (r['bo'] + wor*water_int.loc[i,'bw'])
                gas_cum = oil_cum * r['rs']
                water_cum = wor*oil_cum 
                _df = pd.DataFrame({'np':oil_cum,'gp':gas_cum,'wp':water_cum,'wor':wor,'sw':_sw}, index=[i])
                forecast = forecast.append(_df)

                _so = (1-self.swi)*(1-(oil_cum/self.n))*(r['bo']/ic['bo'])
                _sw = 1-_so
        
        return forecast


                

                


                



