import numpy as np
import pandas as pd 
from .grid import grid
from ...pvtpy.black_oil import oil, gas, water
from ...krpy import water_oil_kr, gas_oil_kr
from ...wellpy.path import wells_group


class numerical:
    def __init__(self, **kwargs):
        self.relaxation = kwargs.pop('relaxation',1)
        self.max_iter = kwargs.pop('max_iter',30)
        self.date_range = kwargs.pop('date_range', None)

    @property
    def relaxation(self):
        return self._relaxation

    @relaxation.setter 
    def relaxation(self,value):
        assert isinstance(value,(int,float))
        assert value >= 0 
        self._relaxation = value 

    @property 
    def max_iter(self):
        return self._max_iter

    @max_iter.setter 
    def max_iter(self,value):
        assert isinstance(value,int) and value > 0
        self._max_iter = value 
    
    @property
    def date_range(self):
        return self._date_range 
    
    @date_range.setter
    def date_range(self,value):
        assert isinstance(value,np.ndarray)
        assert all(isinstance(i,np.datetime64) for i in value)
        self._date_range = value

class sim_model:

    def __init__(self,**kwargs):

        #Grid. Includes petrophysical properties
        self.grid = kwargs.pop('grid',None)

        #Number of phases to simulate
        self.phase = kwargs.pop('phase',None)

        #pvt
        self.pvt = kwargs.pop('pvt',None)

        # kr and pc
        self.rock_fluid = kwargs.pop('rock_fluid', None)

        #Wells
        self.wells = kwargs.pop('wells',None)

        #Numerical
        self.numerical = kwargs.pop('numerical',None)

    ## Properties

    @property
    def grid(self):
        return self._grid

    @grid.setter 
    def grid(self,value):
        assert isinstance(value,grid), f"{type(value)} not allowed"
        assert all( i in list(value.petrophysics.keys()) for i in ['PORO','PERMX','PERMY','PERMZ','RT'])
        self._grid = value

    @property 
    def phase(self):
        return self._phase

    @phase.setter
    def phase(self,value):
        phases = ['oil','water','gas']
        assert isinstance(value,list) and len(value) <= 3
        assert all(i in phases for i in value)
        self._phase = value

    @property 
    def pvt(self):
        return self._pvt 
    
    @pvt.setter
    def pvt(self,value):
        assert isinstance(value,dict)
        
        for i in value:
            assert i in self.phase
        self._pvt  = value

    @property 
    def rock_fluid(self):
        return self._rock_fluid

    @rock_fluid.setter
    def rock_fluid(self,value):
        
        if len(self.phase) == 1:
            self._rock_fluid = None
        else:
            assert isinstance(value,dict)
            rt_list =  np.unique(self.grid.petrophysics['RT']).tolist()

            #Assert all rock types are present in rock fluid
            assert all(i in list(value.keys()) for i in rt_list)

            for i in value:
                #Assert the keys for each rocktype are dictionary
                assert isinstance(value[i], dict)

                if len(self.phase) == 2: 
                    assert len(list(value[i].keys())) == 1
                else: 
                    assert len(list(value[i].keys())) == 2

                for j in value[i]:
                    assert j in ['krwo','krgo']
                    assert isinstance(value[i][j],(water_oil_kr,gas_oil_kr))

            self._rock_fluid = value

    
    @property
    def wells(self):
        return self._wells 
    
    @wells.setter
    def wells(self,value):
        if value is not None:
            assert isinstance(value, wells_group)
        self._wells = value

    @property
    def numerical(self):
        return self._numerical
    
    @numerical.setter
    def numerical(self,value):
        assert isinstance(value, numerical)
        self._numerical = value
    