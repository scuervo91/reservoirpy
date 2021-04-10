import pandas as pd 
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def kr_curve(sn:np.ndarray, n:float, krend:float) -> np.ndarray:
    """kr_curve [Estimate Relative permeability curve from normalized saturation, exponent and end-points]

    Parameters
    ----------
    sn : np.ndarray
        [Saturation Array]
    n : float
        [Exponent]
    krend : float
        [End-point]

    Returns
    -------
    np.ndarray
        [Relative permeability curve]
    """
    return krend * np.power(sn,n)

def sw_normalize(sw:np.ndarray, swir:float, sor:float) -> np.ndarray:
    """sw_normalize [Convert array of water saturation to normalized water saturation]

    Parameters
    ----------
    sw : np.ndarray
        [Water saturation array]
    swir : float
        [Irreducible water saturation]
    sor : float
        [Residual Oil Saturation]

    Returns
    -------
    np.ndarray
        [Normalized water Saturation]
    """
    swn = (sw - swir) / (1 - swir - sor)
    return swn

def sw_denormalize(swn:np.ndarray, swir:float, sor:float) -> np.ndarray:
    """sw_normalize [Convert array of normalized water saturation to water saturation]

    Parameters
    ----------
    sw : np.ndarray
        [normalized Water saturation array]
    swir : float
        [Irreducible water saturation]
    sor : float
        [Residual Oil Saturation]

    Returns
    -------
    np.ndarray
        [water Saturation]
    """
    sw = swn * (1 - swir - sor) + swir 
    return sw


class Kr(pd.DataFrame):
    
    def __init__(self, *args, **kwargs):
        wet_col = kwargs.pop("index", 'sw')
        wet = kwargs.pop('wet',None)
        assert wet_col in ['sw','sl','so']
        assert isinstance(wet,(list,np.ndarray,type(None)))
        super().__init__(*args, **kwargs)
                
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
        return Kr

class KrWaterOil:

    def __init__(self, **kwargs):
        self.swir = kwargs.pop('swir',0)
        self.sor = kwargs.pop('sor',0)
        self.krwend = kwargs.pop('krwend',1)
        self.kroend = kwargs.pop('kroend',1)
        self.pcend = kwargs.pop('pcend',0)
        self.nw = kwargs.pop('nw',1)
        self.no = kwargs.pop('no',1)
        self.np = kwargs.pop('np',1)
        self.kr = kwargs.pop('kr', None)

    #Properties

    @property
    def swir(self):
        return self._swir

    @swir.setter
    def swir(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0 and value <= 1 
        self._swir = value

    @property
    def sor(self):
        return self._sor

    @sor.setter
    def sor(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0 and value <= 1 
        self._sor = value

    @property
    def krwend(self):
        return self._krwend

    @krwend.setter
    def krwend(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0 and value <= 1 
        self._krwend = value

    @property
    def kroend(self):
        return self._kroend

    @kroend.setter
    def kroend(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0 and value <= 1 
        self._kroend = value

    @property
    def pcend(self):
        return self._pcend

    @pcend.setter
    def pcend(self,value):
        assert isinstance(value,(int,float))
        self._pcend = value

    @property
    def nw(self):
        return self._nw

    @nw.setter
    def nw(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0
        self._nw = value

    @property
    def no(self):
        return self._no

    @no.setter
    def no(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0
        self._no = value

    @property
    def np(self):
        return self._np

    @np.setter
    def np(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0
        self._np = value

    @property
    def kr(self):
        return self._kr

    @kr.setter
    def kr(self,value):
        if value is not None:
            assert isinstance(value,Kr)
        self._kr = value

    #Methods

    def build_kr(self, n=10):

        #Make Normalized Sw and So
        swn = np.linspace(0,1,n)
        son = 1 - swn 

        #Calculate Krw, kro  and pc
        kro = kr_curve(son,self.no,self.kroend)
        kro = np.append(kro,0)
        
        krw = kr_curve(swn,self.nw,self.krwend)
        krw = np.append(krw,1)
        
        pcwo = self.pcend * np.power(son,self.np) 
        pcwo = np.append(pcwo,0)
        
        kr_ratio = kro/krw

        #Calculate Sw from endpoints
        sw = sw_denormalize(swn, self.swir, self.sor)
        sw = np.append(sw,1)
        
        col_dict = {
            'sw':sw,
            'krw':krw,
            'kro':kro,
            'pcwo':pcwo,
            'kr_ratio':kr_ratio
        }
        
        kr_table = Kr(col_dict)

        self._kr = kr_table

        return kr_table
    
    def fw(self, muo, bo, muw, bw):
        
        if self.kr is not None:
            kr_table = self.kr
            fw = 1 / (1 + (muw/muo)*self.kr['kr_ratio'])

            kr_table['fw'] = fw 
            self._kr = kr_table
            

    def plot(self,
        ax=None,
        norm=False, 
        ann=False, 
        pc = False,
        kr = True,
        krw_kw={},
        kro_kw={}, 
        ann_kw = {},
        pcwo_kw = {}
        ):

        if self.kr is not None:  

            #Set default plot properties krw
            def_krw_kw = {
                'color': 'blue',
                'linestyle':'--',
                'linewidth': 2
                }    

            for (k,v) in def_krw_kw.items():
                if k not in krw_kw:
                    krw_kw[k]=v

            #Set default plot properties kro
            def_kro_kw = {
                'color': 'green',
                'linestyle':'--',
                'linewidth': 2
                }    

            for (k,v) in def_kro_kw.items():
                if k not in kro_kw:
                    kro_kw[k]=v

            def_pc_kw = {
                'color': 'black',
                'linestyle':'--',
                'linewidth': 1
                }    

            for (k,v) in def_pc_kw.items():
                if k not in pcwo_kw:
                    pcwo_kw[k]=v

            #Set default plot properties kro
            def_ann_kw = {
                'xytext': (0,15),
                'textcoords':'offset points',
                'arrowprops': {'arrowstyle':"->"},
                'bbox':{'boxstyle':'round', 'fc':'0.8'},
                'fontsize':11
                }    

            for (k,v) in def_ann_kw.items():
                if k not in ann_kw:
                    ann_kw[k]=v


            if kr:
                if norm:
                    sw_x = (self.kr.index[:-1] - self.swir) / (1 - self.swir - self.sor)
                    krw = self.kr['krw'].iloc[:-1].values
                    kro = self.kr['kro'].iloc[:-1].values
                else:
                    sw_x = self.kr.index 
                    krw = self.kr['krw'].values
                    kro = self.kr['kro'].values

                #Set the axes      
                krax= ax or plt.gca()
                krax.plot(sw_x, krw, **krw_kw)
                krax.plot(sw_x, kro, **kro_kw)

                #set labels
                krax.set_xlabel('Water Saturation []')
                krax.set_ylabel('Kr []')
                krax.set_xlim([0,1])
                krax.set_ylim([0,1])


            if pc and not norm:
                if kr:
                    pcax=krax.twinx()
                    pcax.yaxis.set_label_position("right")
                    pcax.yaxis.tick_right()
                else:
                    pcax= ax or plt.gca()

                pcax.plot(self.kr.index, self.kr['pcwo'], **pcwo_kw)
                pcax.set_ylabel('Capillary Pressure [psi]')

            #Annotate
            if ann and kr and not norm:
                krax.annotate(
                    'swir',
                    xy = (self.kr.index[0],self.kr['krw'].iloc[0]),
                    xycoords='data',
                    **ann_kw
                ) 
                krax.annotate(
                    'swor',
                    xy = (self.kr.index[-2],self.kr['kro'].iloc[-2]),
                    xycoords='data',
                    **ann_kw
                ) 
                krax.annotate(
                    'kroend',
                    xy = (self.kr.index[0],self.kr['kro'].iloc[0]),
                    xycoords='data',
                    **ann_kw
                ) 
                krax.annotate(
                    'krwend',
                    xy = (self.kr.index[-2],self.kr['krw'].iloc[-2]),
                    xycoords='data',
                    **ann_kw
                ) 
                

        else:
            raise ValueError('kr is not defiend')

    def to_ecl(self):
        
        assert self.kr is not None
        string = "SWOF\n"
        
        string += self.kr.reset_index().to_string(header=False, index=False) +'/\n'
        
        return string

    def fit(self, df:pd.DataFrame, krw:str=None, kro:str=None):
        
        sw = df.index.values 
        
        if krw is not None:
            krw_array = df[krw].values
        
            popt, pcov = curve_fit(kr_curve, sw, krw_array, bounds=([0.01,0], [np.inf, 1]))
            
            print(f'Krw parameters\n-----\n n: {popt[0]}\n krend: {popt[1]}')
            self.nw = popt[0]
            self.krwend = popt[1]
            
        if kro is not None:
            kro_array = df[kro].values
        
            popt, pcov = curve_fit(kr_curve, 1-sw, kro_array, bounds=([0.01,0], [np.inf, 1]))
            
            print(f'Kro parameters\n-----\n n: {popt[0]}\n krend: {popt[1]}')
            self.no = popt[0]
            self.kroend = popt[1]      

    
class KrGasOil:

    def __init__(self, **kwargs):
        self.slc = kwargs.pop('slc',0)
        self.sgc = kwargs.pop('sgc',0)
        self.krgend = kwargs.pop('krgend',1)
        self.kroend = kwargs.pop('kroend',1)
        self.pcend = kwargs.pop('pcend',0)
        self.no = kwargs.pop('no',1)
        self.ng = kwargs.pop('ng',1)
        self.np = kwargs.pop('np',1)
        self.kr = kwargs.pop('kr', None)

    #Properties

    @property
    def slc(self):
        return self._slc

    @slc.setter
    def slc(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0 and value <= 1 
        self._slc = value

    @property
    def sgc(self):
        return self._sgc

    @sgc.setter
    def sgc(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0 and value <= 1 
        self._sgc = value

    @property
    def krgend(self):
        return self._krgend

    @krgend.setter
    def krgend(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0 and value <= 1 
        self._krgend = value

    @property
    def kroend(self):
        return self._kroend

    @kroend.setter
    def kroend(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0 and value <= 1 
        self._kroend = value

    @property
    def pcend(self):
        return self._pcend

    @pcend.setter
    def pcend(self,value):
        assert isinstance(value,(int,float))
        self._pcend = value

    @property
    def ng(self):
        return self._ng

    @ng.setter
    def ng(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0
        self._ng = value

    @property
    def no(self):
        return self._no

    @no.setter
    def no(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0
        self._no = value

    @property
    def np(self):
        return self._np

    @np.setter
    def np(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0
        self._np = value

    @property
    def kr(self):
        return self._kr

    @kr.setter
    def kr(self,value):
        if value is not None:
            assert isinstance(value,Kr)
        self._kr = value

    #Methods

    def build_kr(self, n=10):

        #Make Normalized Sw and So
        sgn = np.linspace(0,1,n)
        son = 1 - sgn 

        #Calculate Krw, kro  and pc
        kro = self.kroend * np.power(son,self.no)
        krg = self.krgend * np.power(sgn,self.ng)
        pcgo = self.pcend * np.power(sgn,self.np) 

        #Calculate Sg from endpoints
        sg = sgn * (1 - self.slc - self.sgc) + self.sgc
        sl = 1-sg
        kr_table = Kr({
            'sl':sl,
            'sg':sg,
            'krg':krg,
            'kro':kro,
            'pcgo':pcgo,
        }, index='sg')

        self._kr = kr_table

        return kr_table

    def plot(self,
        ax=None,
        ann=False, 
        pc = False,
        kr = True,
        krg_kw={},
        kro_kw={}, 
        ann_kw = {},
        pcgo_kw = {}
        ):

        if self.kr is not None:  

            #Set default plot properties krw
            def_krg_kw = {
                'color': 'red',
                'linestyle':'--',
                'linewidth': 2
                }    

            for (k,v) in def_krg_kw.items():
                if k not in krg_kw:
                    krg_kw[k]=v

            #Set default plot properties kro
            def_kro_kw = {
                'color': 'gray',
                'linestyle':'--',
                'linewidth': 2
                }    

            for (k,v) in def_kro_kw.items():
                if k not in kro_kw:
                    kro_kw[k]=v

            def_pc_kw = {
                'color': 'black',
                'linestyle':'--',
                'linewidth': 1
                }    

            for (k,v) in def_pc_kw.items():
                if k not in pcgo_kw:
                    pcgo_kw[k]=v

            #Set default plot properties kro
            def_ann_kw = {
                'xytext': (0,15),
                'textcoords':'offset points',
                'arrowprops': {'arrowstyle':"->"},
                'bbox':{'boxstyle':'round', 'fc':'0.8'},
                'fontsize':11
                }    

            for (k,v) in def_ann_kw.items():
                if k not in ann_kw:
                    ann_kw[k]=v


            if kr:
                sw_x = self.kr.index 
                krg = self.kr['krg'].values
                kro = self.kr['kro'].values

                #Set the axes      
                krax= ax or plt.gca()
                krax.plot(sw_x, krg, **krg_kw)
                krax.plot(sw_x, kro, **kro_kw)

                #set labels
                krax.set_xlabel('Gas Saturation []')
                krax.set_ylabel('Kr []')
                krax.set_xlim([0,1])
                krax.set_ylim([0,1])


            if pc:
                if kr:
                    pcax=krax.twinx()
                    pcax.yaxis.set_label_position("right")
                    pcax.yaxis.tick_right()
                else:
                    pcax= ax or plt.gca()

                pcax.plot(self.kr.index, self.kr['pcwo'], **pcgo_kw)
                pcax.set_ylabel('Capillary Pressure [psi]')

            #Annotate
            if ann and kr:
                krax.annotate(
                    'slc',
                    xy = (self.kr.index[0],self.kr['kro'].iloc[0]),
                    xycoords='data',
                    **ann_kw
                ) 
                krax.annotate(
                    'sgc',
                    xy = (self.kr.index[-1],self.kr['krg'].iloc[-1]),
                    xycoords='data',
                    **ann_kw
                ) 
                krax.annotate(
                    'kroend',
                    xy = (self.kr.index[0],self.kr['kro'].iloc[0]),
                    xycoords='data',
                    **ann_kw
                ) 
                krax.annotate(
                    'krgend',
                    xy = (self.kr.index[-1],self.kr['krw'].iloc[-1]),
                    xycoords='data',
                    **ann_kw
                ) 
                
        else:
            raise ValueError('kr is not defiend')
        
    def to_ecl(self):
        
        assert self.kr is not None
        string = "SGOF\n"
        
        string += self.kr.reset_index().to_string(header=False, index=False) +'/\n'
        
        return string

    def fit(self, df:pd.DataFrame, krg:str=None, kro:str=None):
        
        sg = df.index.values 
        
        if krg is not None:
            krg_array = df[krg].values
        
            popt, pcov = curve_fit(kr_curve, sg, krg_array, bounds=([0.01,0], [np.inf, 1]))
            
            print(f'Krg parameters\n-----\n n: {popt[0]}\n krend: {popt[1]}')
            self.nw = popt[0]
            self.krgend = popt[1]
            
        if kro is not None:
            kro_array = df[kro].values
        
            popt, pcov = curve_fit(kr_curve, 1-sg, kro_array, bounds=([0.01,0], [np.inf, 1]))
            
            print(f'Kro parameters\n-----\n n: {popt[0]}\n krend: {popt[1]}')
            self.no = popt[0]
            self.kroend = popt[1]     