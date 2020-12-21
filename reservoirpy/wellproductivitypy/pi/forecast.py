import numpy as np
import pandas as pd 
from ...wellpy.path import Well, WellsGroup
from .inflow import OilInflow, GasInflow
from ...pvtpy.black_oil import Oil, Gas, Water


class pressure_model:
    def __init__(self,**kwargs):

        self.b = kwargs.pop('b',0)
        self.x_oil = kwargs.pop('x_oil',0)
        self.x_gas = kwargs.pop('x_gas',0)

    def get_value(self,oil_cum=0, gas_cum=0):
        x = np.array([1,oil_cum,gas_cum])
        theta = np.array([self.b,self.x_oil,self.x_gas])
        p = np.dot(x,theta)
        return p


class forecast_model:

    def __init__(self, **kwargs):

        self.wells_group = kwargs.pop('wells_group', None)
        self.formations = kwargs.pop('formations',None)
        self.initial_conditions = kwargs.pop('initial_conditions',None)
        self.pressure_model = kwargs.pop('pressure_model',None)
        self.wells_restrictions = kwargs.pop('wells_restrictions',None)
        self.group_restrictions = kwargs.pop('group_restrictions',None)
        self.fluids = kwargs.pop('fluids',None)
        self.productivity_index = kwargs.pop('productivity_index',None)

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

    @property 
    def pressure_model(self):
        return self._pressure_model

    @pressure_model.setter
    def pressure_model(self,value):
        assert isinstance(value,pressure_model)
        self._pressure_model = value 

    @property 
    def productivity_index(self):
        return self._productivity_index

    @productivity_index.setter
    def productivity_index(self,value):
        assert isinstance(value,dict)
        self._productivity_index = value 