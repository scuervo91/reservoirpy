import numpy as np 
import pandas as pd 
from datetime import date, timedelta


def bsw_to_wor(bsw):
    assert isinstance(bsw,(int,float,np.ndarray))
    bsw = np.atleast_1d(bsw)
    assert np.all((bsw>=0)&(bsw<=1))
    wor = bsw/(1-bsw)
    return wor 

def wor_to_bsw(wor):
    assert isinstance(wor,(int,float,np.ndarray))
    wor = np.atleast_1d(wor)
    assert np.all(wor>=0)
    bsw = wor/(wor+1)
    return bsw   


def wor_forecast(
    range_time,
    fluid_rate, 
    slope, wor_i, 
    econ_limit = None,
    np_limit=None, 
    wor_limit=None
    ):
    """
    Estimate a Forecast curve given wor+1 parameters

    Attributes:

        Return -> 
    """
    assert isinstance(range_time,pd.Series)
    days_number = range_time.apply(lambda x: x.toordinal()).values
    assert isinstance(fluid_rate,(pd.Series,np.ndarray))
    fluid_rate = np.atleast_1d(fluid_rate)

    assert fluid_rate.shape == days_number.shape, 'range time and fluid must be the same shape'
    assert isinstance(slope,(int,float))
    assert isinstance(wor_i,(int,float)) and wor_i >= 0

    wor_i1 = wor_i + 1

    df = pd.DataFrame()
    
    cum = 0
    days_diff = np.diff(days_number,prepend=days_number[0]-1)

    for i in range(days_number.shape[0]):
        wor_1 = np.exp(slope*cum)*wor_i1
        wor = wor_1 - 1
        bsw = wor_to_bsw(wor)
        qo = fluid_rate[i]*(1-bsw)
        qw = fluid_rate[i]*bsw
        cum += qo*days_diff[i]
        _df = pd.DataFrame({'qf':fluid_rate[i],'qo':qo,'qw':qw,'bsw':bsw,'wor_1':wor_1,'wor':wor,'np':cum})
        df = df.append(_df)

        if econ_limit is not None:
            if qo <= econ_limit:
                break
        
        if np_limit is not None:
            if cum >= np_limit:
                break     

        if wor_limit is not None:
            if wor >= wor_limit:
                break  

    df.index = range_time[:df.shape[0]]
    df.index.name = 'time'

    return df

class wor_declination:
    def __init__(self,**kwargs):
        self.slope = kwargs.pop('slope',0)
        self.bsw_i = kwargs.pop('bsw_i',0)
        self.start_date = kwargs.pop('start_date',None)
        self.end_date = kwargs.pop('end_date',None)
        self.econ_limit = kwargs.pop('econ_limit', None)
        self.wor_limit = kwargs.pop('wor_limit', None)
        self.np_limit = kwargs.pop('np_limit', None)
        self.fluid_rate = kwargs.pop('fluid_rate', None)

    ## Properties
    @property
    def slope(self):
        return self._slope
    
    @slope.setter
    def slope(self, value):
        assert isinstance(value,(int, float)), 'slope must be a numerical'
        self._slope = value

    @property
    def bsw_i(self):
        return self._bsw_i

    @bsw_i.setter
    def bsw_i(self,value):
        assert isinstance(value,(int, float)), 'bsw must be a numerical'
        assert value >= 0 and value <= 1, 'bsw must be between 0 and 1'
        self._bsw_i = value
        self._wor_i = value/(1-value)

    @property
    def wor_i(self):
        return self._wor_i

    @property
    def start_date(self):
        return self._start_date
    
    @start_date.setter
    def start_date(self,value):
        if value is not None:
            assert isinstance(value,date), f'{type(value)} not accepted. It must be date'
        self._start_date = value

    @property
    def end_date(self):
        return self._end_date

    @end_date.setter
    def end_date(self,value):
        if value is not None:
            assert isinstance(value,date), f'{type(value)} not accepted. It must be date'
        self._end_date = value
  
    @property
    def econ_limit(self):
        return self._econ_limit

    @econ_limit.setter
    def econ_limit(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be number'
        self._econ_limit = value

    @property
    def np_limit(self):
        return self._np_limit

    @np_limit.setter
    def np_limit(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be number'
        self._np_limit = value

    @property
    def wor_limit(self):
        return self._wor_limit

    @wor_limit.setter
    def wor_limit(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be number'
        self._wor_limit = value

    @property
    def fluid_rate(self):
        return self._fluid_rate

    @fluid_rate.setter
    def fluid_rate(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be number'
        self._fluid_rate = value


    def forecast(self,
        time_range:pd.Series=None,
        start_date:date=None, 
        end_date:date=None, 
        fluid_rate:(pd.Series,np.ndarray,int,float)=None, 
        fq:str='M', 
        np_limit=None,
        wor_limit=None,
        econ_limit=None,
        **kwargs
    ):
        """
        Forecast curve from the declination object. 
    
        Input: 

        Return: 

        """
        if time_range is not None:
            assert isinstance(time_range,(pd.Series)), 'start_date must be pd.Series with dates'
        else:
            if start_date is None: 
                assert self.start_date is not None
                start_date = self.start_date
            else:
                assert isinstance(start_date,date), 'start_date must be date'

            if end_date is None: 
                assert self.end_date is not None
                end_date = self.end_date
            else:
                assert isinstance(end_date,date), 'send_date must be date'

            time_range = pd.Series(pd.date_range(start=start_date, end=end_date, freq=fq, **kwargs))

       
        if fluid_rate is not None:
            assert isinstance(fluid_rate,(pd.Series,np.ndarray,int,float))
            fluid_rate = np.atleast_1d(fluid_rate)

            if fluid_rate.shape==(1,):
                fluid_rate = np.full(time_range.shape,fluid_rate)
            else:
                assert fluid_rate.shape == time_range.shape
        else:
            fluid_rate = np.full(time_range.shape,self.fluid_rate)


        if econ_limit is None:
            econ_limit = self.econ_limit

        if np_limit is None:
            np_limit = self.np_limit

        if wor_limit is None:
            wor_limit = self.wor_limit

        f = wor_forecast(
            time_range,
            fluid_rate, 
            self.slope, 
            self.wor_i, 
            np_limit = np_limit, 
            wor_limit = wor_limit, 
            econ_limit = econ_limit
            )
        
        Np = f['np'].iloc[-1]

        return f, Np