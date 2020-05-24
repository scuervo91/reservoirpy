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


def wor_forecast(range_time,fluid_rate, slope, wor_i):
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

    df.index = range_time

    return df

class wor_declination:
    def __init__(self,**kwargs):
        self.slope = kwargs.pop('slope',0)
        self.bsw_i = kwargs.pop('bsw_i',0)

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

    def forecast(self,time_range=None,start_date=None, end_date=None, fluid_rate=None, fq='M', **kwargs):
        """
        Forecast curve from the declination object. 
    
        Input: 

        Return: 

        """
        assert isinstance(time_range,(pd.Series,type(None))), 'start_date must be pd.Series with dates'
        assert isinstance(start_date,(date,type(None))), 'start_date must be date'
        assert isinstance(end_date,(date,type(None))), 'send_date must be date'
        assert isinstance(fluid_rate,(pd.Series,np.ndarray,int,float))
        fluid_rate = np.atleast_1d(fluid_rate)

        # Create the time range
        if time_range is None:
            if end_date is None: 
                end_date = start_date + timedelta(days=365)
            time_range = pd.Series(pd.date_range(start=start_date, end=end_date, freq=fq, **kwargs))

        if fluid_rate.shape==(1,):
            fluid_rate = np.full(time_range.shape,fluid_rate)
        else:
            assert fluid_rate.shape == time_range.shape

        f = wor_forecast(time_range,fluid_rate, self.slope, self.wor_i)

        return f 