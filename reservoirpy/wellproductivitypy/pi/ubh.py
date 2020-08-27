import numpy as np 
import pandas as pd 
from .outflow import incompressible_pressure_profile, two_phase_pressure_profile, two_phase_outflow_curve
from .als import als
from .inflow import oil_inflow



class ubh(als):
    def __init__(self,**kwargs):
        self.n = kwargs.pop('n',15)
        self.depth = kwargs.pop('depth',None)
        self.reference = kwargs.pop('ref', 'tvd')
        self.brand = kwargs.pop('brand',None)
        self.nozzle = kwargs.pop('nozzle',None)
        self.throat = kwargs.pop('throat', None)
        self.power_fluid_ge =kwargs.pop('power_fluid_ge',0.433)
        self.injection_di = kwargs.pop('injection_di',2.99)
        self.return_di = kwargs.pop('return_di',5)
        super().__init__(**kwargs)

    @property 
    def depth(self):
        return self._depth

    @depth.setter 
    def depth(self,value):
        if isinstance(value,(int,float)):
            value = np.full(self.n,value)
        elif isinstance(value,(list,np.ndarray)):
            value = np.atleast_1d(value)
            assert value.ndim == 1 
        else:
            raise ValueError('Type not allowed')
        self._depth = value

    @property
    def reference(self):
        return self._reference 

    @reference.setter 
    def reference(self,value):
        assert isinstance(value,str) and value in ['md', 'tvd']
        self._reference = value

    @property
    def brand(self):
        return self._brand 

    @brand.setter 
    def brand(self, value):
        """brand [Ubh brand name]

        :param value: [Name of the UBH brand]
        :type value: [str]
        """
        if value is not None:
            assert isinstance(value,str)
        self._brand = value 

    @property
    def nozzle(self):
        return self._nozzle 

    @nozzle.setter 
    def nozzle(self,value):
        """nozzle [Nozzle Size mm]

        :param value: [nozzle size mm]
        :type value: [float]
        """

        assert isinstance(value,(int,float))
        assert value > 0
        self._nozzle = value

    @property
    def throat(self):
        return self._throat 

    @throat.setter 
    def throat(self,value):
        """throat [Throat size mm]

        :param value: [Throat size mm]
        :type value: [float]
        """
        assert isinstance(value,(int,float))
        assert value > 0
        self._throat = value

    @property
    def power_fluid_ge(self):
        return self._power_fluid_ge 

    @power_fluid_ge.setter 
    def power_fluid_ge(self,value):
        """power_fluid_ge [power_fluid_ge]

        :param value: [power_fluid_ge ]
        :type value: [float]
        """
        assert isinstance(value,(int,float))
        assert value > 0
        self._power_fluid_ge = value

    @property
    def injection_di(self):
        return self._injection_di 

    @injection_di.setter 
    def injection_di(self,value):
        """injection_di [injection_di]

        :param value: [injection_di ]
        :type value: [float]
        """
        if isinstance(value,(int,float)):
            value = np.full(self.n,value)
        elif isinstance(value,(list,np.ndarray)):
            value = np.atleast_1d(value)
        self._injection_di = value

    @property
    def return_di(self):
        return self._return_di 

    @return_di.setter 
    def return_di(self,value):
        """return_di [return_di]

        :param value: [return_di ]
        :type value: [float]
        """
        if isinstance(value,(int,float)):
            value = np.full(self.n,value)
        elif isinstance(value,(list,np.ndarray)):
            value = np.atleast_1d(value)
        self._return_di = value

    # Methods
    def get_area(self,value):
        assert isinstance(value,str) and value in ['nozzle', 'throat']

        size = self.nozzle if value == 'nozzle' else self.throat

        size_in = size / 25.4 
        area_in_sq = np.pi * np.power(size_in/2,2)

        return area_in_sq

    def fad(self):
        nozzle_area = self.get_area('nozzle')
        throat_area = self.get_area('throat')

        fad = nozzle_area / throat_area 

        return fad
    
    def annulus_area(self):
        nozzle_area = self.get_area('nozzle')
        throat_area = self.get_area('throat')

        annulus_area = throat_area - nozzle_area 

        return annulus_area
"""
    def flow_match(
        self,
        thp=None,
        oil_obj=None,
        gas_obj=None,
        water_obj=None,



    )

"""