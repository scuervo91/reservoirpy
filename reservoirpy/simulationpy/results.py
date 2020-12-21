import numpy as np

class Results:
    def __init__(self,**kwargs):
        
        self.time = kwargs.pop('time',None)
        self.ooip = kwargs.pop('ooip', None)
        self.ogip = kwargs.pop('ogip', None)
        self.owip = kwargs.pop('owip', None)
        self.pore_volume = kwargs.pop('pore_volume', None)
        self.pw = kwargs.pop('pw', None)
        self.po = kwargs.pop('po', None)
        self.pg = kwargs.pop('pg', None)
        self.sw = kwargs.pop('sw', None)
        self.so = kwargs.pop('so', None)
        self.sg = kwargs.pop('sg', None)
        self.qw = kwargs.pop('qw', None)
        self.qo = kwargs.pop('qo', None)
        self.qg = kwargs.pop('qg', None)
        self.bhp = kwargs.pop('bhp', None)

    #Properties

    @property
    def time(self):
        return self._time
    
    @time.setter
    def time(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 1
        self._time = value

    @property
    def ooip(self):
        return self._ooip
    
    @ooip.setter
    def ooip(self,value):
        assert isinstance(value,(int, float, np.int64, np.float64))
        self._ooip = value

    @property
    def ogip(self):
        return self._ogip
    
    @ogip.setter
    def ogip(self,value):
        assert isinstance(value,(int, float, np.int64, np.float64))
        self._ogip = value

    @property
    def owip(self):
        return self._owip
    
    @owip.setter
    def owip(self,value):
        assert isinstance(value,(int, float, np.int64, np.float64))
        self._owip = value

    @property
    def pore_volume(self):
        return self._pore_volume
    
    @pore_volume.setter
    def pore_volume(self,value):
        assert isinstance(value,(int, float, np.int64, np.float64))
        self._pore_volume = value

    @property
    def pw(self):
        return self._pw
    
    @pw.setter
    def pw(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 2
        self._pw = value

    @property
    def po(self):
        return self._po
    
    @po.setter
    def po(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 2
        self._po = value

    @property
    def pg(self):
        return self._pg
    
    @pg.setter
    def pg(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 2
        self._pg = value

    @property
    def sw(self):
        return self._sw
    
    @sw.setter
    def sw(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 2
        self._sw = value

    @property
    def so(self):
        return self._so
    
    @so.setter
    def so(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 2
        self._so = value

    @property
    def sg(self):
        return self._sg
    
    @sg.setter
    def sg(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 2
        self._sg = value

    @property
    def qw(self):
        return self._qw
    
    @qw.setter
    def qw(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 2
        self._qw = value

    @property
    def qo(self):
        return self._qo
    
    @qo.setter
    def qo(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 2
        self._qo = value

    @property
    def qg(self):
        return self._qg
    
    @qg.setter
    def qg(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 2
        self._qg = value

    @property
    def bhp(self):
        return self._bhp
    
    @bhp.setter
    def bhp(self,value):
        assert isinstance(value,np.ndarray) and value.ndim == 2
        self._bhp = value