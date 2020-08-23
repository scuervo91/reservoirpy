import numpy as np

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