import numpy as np 
import pyvista as pv 
import vtk
from shapely.geometry import Point
import math

petrophysical_properties = ['PORO','PERMX','PERMY','PERMZ','SW']

## Auxilary functions
def cell_id(i,j,k,nx,ny):
    cell = (nx*(j-1)+i)+(k-1)*nx*ny
    return cell

def cell_ijk(cell_id,nx,ny):
    k=math.ceil(cell_id/(nx*ny))
    j=math.ceil((cell_id-(nx*ny)*(k-1))/nx)
    i=math.ceil(cell_id-(nx*ny*(k-1))-nx*(j-1))
    return i,j,k


## Grid Class



class grid():
    """
    Class for Reservoir Simulation Grid 
        * Cartesian
        * Corner-Point
    """
    
    def __init__(self,**kwargs):
        #Grid Type
        self._grid_type = kwargs.pop('grid_type',None)

        # Dimensions
        self.nx = kwargs.pop('nx',None)
        self.ny = kwargs.pop('ny',None)
        self.nz = kwargs.pop('nz',None)

        #Cartesian Grid
        self.dx = kwargs.pop('dx',None)
        self.dy = kwargs.pop('dy',None)
        self.dz = kwargs.pop('dz',None)
        self.tops = kwargs.pop('tops',None)
        self.origin = kwargs.pop('origin',None)

        #Corner Point Grid
        self.coord = kwargs.pop('coord',None)
        self.zcorn = kwargs.pop('zcorn',None)

        #Petrophysical properties
        self._petrophysics = {}
        self.petrophysics = kwargs.pop('petrophysics',None)

#####################################################
############## Properties ###########################

    #---Grid Type---
    @property
    def grid_type(self):
        return self._grid_type

    @grid_type.setter
    def grid_type(self,value):
        assert isinstance(value,str), f'{type(value)} not accepted. Name must be str'
        assert value in ['cartesian','corner_point']
        self._grid_type = value

    #---NX---
    @property
    def nx(self):
        return self._nx

    @nx.setter 
    def nx(self,value):
        assert isinstance(value,int) and value > 0 
        self._nx = value 

    #---NY---
    @property
    def ny(self):
        return self._ny

    @ny.setter 
    def ny(self,value):
        assert isinstance(value,int) and value > 0 
        self._ny = value 

    #---NZ---
    @property
    def nz(self):
        return self._nz

    @nz.setter 
    def nz(self,value):
        assert isinstance(value,int) and value > 0 
        self._nz = value 

    #---N---
    @property
    def n(self):
        return self._nx * self._ny * self._nz

    # Cartesian Grid

    #---DX---
    @property
    def dx(self):
        assert len(self._dx) == self.n
        return self._dx 

    @dx.setter 
    def dx(self,value):
        assert isinstance(value,(int,float,list,np.ndarray,type(None))), 'Value entered is not the right type Tyepes allowed: (int,float,list,np.ndarray,None)'
        
        if value is None:
            assert self._grid_type != 'cartesian', 'If cartesian grid set, dx, dy, dz must be set'
            self._dx = value
        elif isinstance(value, (int,float)) and value > 0:
            self._dx = np.full(self.n,value)
        elif isinstance(value,list):
            _dx = np.array(value).flatten(order='F')
            assert len(_dx) == self.n, f'list must be of length {self.n}'
            assert np.issubdtype(_dx.dtype, np.number), 'List must contain only numbers'
            assert all(_dx>0), 'Deltas must be greater than 0'
            self._dx = _dx
        elif isinstance(value,np.ndarray):
            _dx = value.flatten(order='F')
            assert len(_dx) == self.n, f'list must be of length {self.n}'
            assert np.issubdtype(_dx.dtype, np.number), 'List must contain only numbers'
            assert all(_dx>0), 'Deltas must be greater than 0'
            self._dx = _dx

    #---DY---
    @property
    def dy(self):
        assert len(self._dx) == self.n
        return self._dy 

    @dy.setter 
    def dy(self,value):
        assert isinstance(value,(int,float,list,np.ndarray,type(None))), 'Value entered is not the right type Tyepes allowed: (int,float,list,np.ndarray,None)'
        
        if value is None:
            assert self._grid_type != 'cartesian', 'If cartesian grid set, dx, dy, dz must be set'
            self._dy = value
        elif isinstance(value, (int,float)) and value > 0:
            self._dy = np.full(self.n,value)
        elif isinstance(value,list):
            _dy = np.array(value).flatten(order='F')
            assert len(_dy) == self.n, f'list must be of length {self.n}'
            assert np.issubdtype(_dy.dtype, np.number), 'List must contain only numbers'
            assert all(_dy>0), 'Deltas must be greater than 0'
            self._dy = _dy
        elif isinstance(value,np.ndarray):
            _dy = value.flatten(order='F')
            assert len(_dy) == self.n, f'list must be of length {self.n}'
            assert np.issubdtype(_dy.dtype, np.number), 'List must contain only numbers'
            assert all(_dy>0), 'Deltas must be greater than 0'
            self._dy = _dy

    #---DZ---
    @property
    def dz(self):
        assert len(self._dx) == self.n
        return self._dz 

    @dz.setter 
    def dz(self,value):
        assert isinstance(value,(int,float,list,np.ndarray,type(None))), 'Value entered is not the right type Tyepes allowed: (int,float,list,np.ndarray,None)'
        
        if value is None:
            assert self._grid_type != 'cartesian', 'If cartesian grid set, dx, dy, dz must be set'
            self._dz = value
        elif isinstance(value, (int,float)) and value > 0:
            self._dz = np.full(self.n,value)
        elif isinstance(value,list):
            _dz = np.array(value).flatten(order='F')
            assert len(_dz) == self.n, f'list must be of length {self.n}'
            assert np.issubdtype(_dz.dtype, np.number), 'List must contain only numbers'
            assert all(_dz>0), 'Deltas must be greater than 0'
            self._dz = _dz
        elif isinstance(value,np.ndarray):
            _dz = value.flatten(order='F')
            assert len(_dz) == self.n, f'list must be of length {self.n}'
            assert np.issubdtype(_dz.dtype, np.number), 'List must contain only numbers'
            assert all(_dz>0), 'Deltas must be greater than 0'
            self._dz = _dz

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self,value):
        assert isinstance(value,(Point,type(None))), 'Origin point must be Point or None types'
        if value is None:
            assert self._grid_type != 'cartesian', 'If cartesian grid set, origin must be set'
        else:
            assert value.has_z, 'Point must have x,y,z coordinates'
        self._origin = value

    # Corner Point Grid
    #---COORD---
    @property
    def coord(self):
        return self._coord 

    @coord.setter 
    def coord(self,value):
        assert isinstance(value,(list,np.ndarray,type(None))), 'Origin point must be Point or None types'
   
        if value is None:
            assert self.grid_type != 'corner_point', 'If Corner Point Grid set, coord must be set'
            self._coord = value 
        elif isinstance(value,list):
            _coord = np.array(value).flatten(order='F')
            assert len(_coord) == 6*(self.nx+1)*(self.ny+1), f'list must be of length {6*(self.nx+1)*(self.ny+1)}'
            assert np.issubdtype(_coord.dtype, np.number), 'List must contain only numbers'
            self._coord = _coord
        elif isinstance(value,np.ndarray):
            _coord = value.flatten(order='F')
            assert len(_coord) == 6*(self.nx+1)*(self.ny+1), f'list must be of length {6*(self.nx+1)*(self.ny+1)}'
            assert np.issubdtype(_coord.dtype, np.number), 'List must contain only numbers'
            self._coord = _coord

    #---ZCORN---
    @property
    def zcorn(self):
        return self._zcorn 

    @zcorn.setter 
    def zcorn(self,value):
        assert isinstance(value,(list,np.ndarray,type(None))), 'Origin point must be Point or None types'
   
        if value is None:
            assert self.grid_type != 'corner_point', 'If Corner Point Grid set, zcorn must be set'
            self._zcorn = value 
        elif isinstance(value,list):
            _zcorn = np.array(value).flatten(order='F')
            assert len(_zcorn) == 8*self.n, f'list must be of length {8*self.n}'
            assert np.issubdtype(_zcorn.dtype, np.number), 'List must contain only numbers'
            self._zcorn = _zcorn
        elif isinstance(value,np.ndarray):
            _zcorn = value.flatten(order='F')
            assert len(_zcorn) == 8*self.n, f'list must be of length {8*self.n}'
            assert np.issubdtype(_zcorn.dtype, np.number), 'List must contain only numbers'
            self._zcorn = _zcorn

    @property 
    def petrophysics(self):
        return self._petrophysics 

    @petrophysics.setter 
    def petrophysics(self,value):
        assert isinstance(value,(dict,type(None)))

        if isinstance(value,dict):
            for i in value:
                assert i in petrophysical_properties, f"Keyword {i} not in supported properties {petrophysical_properties} "
                assert isinstance(value[i],(int,float,list,np.ndarray))     

                if isinstance(value[i],(int,float)):
                    self._petrophysics[i] = np.full(self.n,value[i])
                elif isinstance(value[i],list):
                    _prop = np.array(value[i]).flatten(order='F')
                    assert len(_prop) == self.n, f'{i} list must be of length {self.n}'
                    assert np.issubdtype(_prop.dtype, np.number), f'{i} List must contain only numbers'
                    assert all(_prop>=0), f'{i} must be greater than 0'
                    self._petrophysics[i] = _prop
                elif isinstance(value[i],np.ndarray):
                    _prop = value[i].flatten(order='F')
                    assert len(_prop) == self.n, f'{i} list must be of length {self.n}'
                    assert np.issubdtype(_prop.dtype, np.number), f'{i} List must contain only numbers'
                    assert all(_prop>=0), f'{i} must be greater than 0'
                    self._petrophysics[i] = _prop

#####################################################
############## Methods ###########################


                






    
    

    

    




