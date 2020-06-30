import numpy as np
import pandas as pd 
from ...wellpy.path import well, wells_group
from inflow import gas_inflow, oil_inflow
from ...pvtpy.black_oil import oil, gas, water


class pressure_model:
    def __init__(self,**kwargs):

        self.b = kwargs.pop('b',0)
        self.oil_cum = kwargs.pop('oil_cum',0)
        self.gas_cum = kwargs.pop('gas_cum',0)

    @property
    def b(self):
        return self._b

    @b.setter
    def b(self,value):
        assert isinstance(value,(int,float))
        self._b = value

    @property
    def oil_cum(self):
        return self._oil_cum

    @oil_cum.setter
    def oil_cum(self,value):
        assert isinstance(value,(int,float))
        self._oil_cum = value

    @property
    def gas_cum(self):
        return self._gas_cum

    @gas_cum.setter
    def gas_cum(self,value):
        assert isinstance(value,(int,float))
        self._gas_cum = value

    def get_value(self,oil_cum=0, gas_cum=0):
        x = np.array([1,oil_cum,gas_cum])
        theta = np.array([self.b,self.oil_cum,self.gas_cum])
        p = np.dot(x,theta)
        return p

class pi_forecast_model:

    def __init__(self, **kwargs):
        
        self.wells_group = kwargs.pop('wells_group', None)
        self.formations = kwargs.pop('formations',None)
        self.initial_conditions = kwargs.pop('initial_conditions',None)
        self.params = kwargs.pop('params',None)
        self.wells_restrictions = kwargs.pop('wells_restrictions',None)
        self.group_restrictions = kwargs.pop('group_restrictions',None)
        self.fluids = kwargs.pop('fluids',None)

    @property
    def wells_group(self):
        return self._wells_group

    @wells_group.setter
    def wells_group(self, value):
        assert isinstance(value, wells_group)
        self._wells_group = value

    @property
    def formations(self):
        return self._formations
    
    @formations.setter
    def formations(self,value):
        assert isinstance(value,(str,list))
        if isinstance(value,str):
            value = [value]
        self._formations = value

    @property 
    def params(self):
        return self._params

    @params.setter
    def params(self,value):
        assert isinstance(value,pressure_model)
        self._params = value 
    
    @property
    def initial_conditions(self):
        return self._initial_conditions

    @initial_conditions.setter
    def initial_conditions(self,value):
        assert isinstance(value,dict)
        assert len(self.formations) == len(value)

        _keys = ['oil_cum','gas_cum','water_cum','pi']

        for i in value:
            assert i in self.formations
            assert isinstance(value[i],dict) 
            for k in value[i]:
                assert k in _keys
                assert isinstance(value[i][k],(int,float))

        self._initial_conditions = value

    @property
    def fluids(self):
        return self._fluids

    @fluids.setter
    def fluids(self,value):
        assert isinstance(value,(dict,type(None)))
        if isinstance(value,dict):
            assert len(self.formations) == len(value)
            for i in value:
                assert i in self.formations
                assert isinstance(value[i],dict) 
                for k in value[i]:
                    assert isinstance(value[i][k],(oil,gas,water))

        self._fluids = value



    


    

