import numpy as np
import pandas as pd 
from .grid import grid
from ...pvtpy.black_oil import oil, gas, water


class sim_model:

    def __init__(self,**kwargs):

        #Grid. Includes petrophysical properties
        self.grid = kwargs.pop('grid',None)

        #Number of phases to simulate
        self.phases = kwargs.pop('phases',None)

        #pvt
        self.oil_model = kwargs.pop('oil_model',None)
        self.water_model = kwargs.pop('water_model',None)
        self.gas_model = kwargs.pop('gas_model',None)

    ## Properties

    @property
    def grid(self):
        return self._grid

    @grid.setter 
    def grid(self,value):
        assert isinstance(value,grid), f"{type(value)} not allowed"
        self._grid = value