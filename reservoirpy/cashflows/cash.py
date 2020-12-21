
import numpy as np
import pandas as pd
from datetime import date, timedelta
from cashflows.timeseries import cashflow


class CashFlow:
    def __init__(self, **kwargs):
        self.const_value = kwargs.pop('const_value',0)
        self.start = kwargs.pop('start',None)
        self.end = kwargs.pop('end',None)
        self.periods = kwargs.pop('end',None)
        self.freq = kwargs.pop('freq','M')
        self.chgpts = kwargs.pop('chgpts',None)
        self.name = kwargs.pop('name', None)

    @property
    def const_value(self):
        return self._const_value

    @const_value.setter
    def const_value(self,value):
        assert isinstance(value,(int,float, list))
        self._const_value = value 

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self,value):
        if value is not None:
            assert isinstance(value,(date,str)), f'{type(value)} not accepted'
        self._start = value 

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self,value):
        if value is not None:
            assert isinstance(value,(date,str))
        self._end = value 

    @property
    def freq(self):
        return self._freq

    @freq.setter
    def freq(self,value):
        assert isinstance(value,str), f"{type(value)} not accepted. Name must be str"
        assert value in ['A', 'BA', 'Q', 'BQ', 'M', 'BM', 'CBM', 'SM', '6M', '6BM', '6CMB']      
        self._freq = value

    @property
    def period(self):
        return self._period

    @period.setter
    def period(self,value):
        if value is not None:
            assert isinstance(value,int)
        self._period = value 

    @property
    def chgpts(self):
        return self._chgpts

    @chgpts.setter
    def chgpts(self,value):
        if value is not None:
            assert isinstance(value,dict)
        self._chgpts = value 

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,value):
        if value is not None:
            assert isinstance(value,str)
        self._name = value 


    def cashflow(self):

        ch = cashflow(
            const_value=self.const_value,
            start = self.start,
            end=self.end, 
            periods=self.periods,
            freq=self.freq,
            chgpts=self.chgpts, 
            name = self.name
        )
        return ch
