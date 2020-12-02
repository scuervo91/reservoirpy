#########################################################################
# Most of the Code has been taken from the next  Github Repository:      #
#  https://github.com/BinWang0213/PyGRDECL                              #
#  Code is used to load and manipulate Eclipse Data Grid
#########################################################################

import numpy as np 
import pyvista as pv 
import vtk
from shapely.geometry import Point
import math
import os
import pandas as pd 

petrophysical_properties = ['PORO','PERMX','PERMY','PERMZ','SW','RT']

SupportKeyWords=[
    'SPECGRID', #Dimenion of the corner point grid
    'DIMENS',   #Define the dimension of the cartesian grid
    'TOPS','DX','DY','DZ',
    'COORD','ZCORN',
    'PORO',
    'PERMX' , 'PERMXY', 'PERMXZ', 
    'PERMYX', 'PERMY' , 'PERMYZ', 
    'PERMZX', 'PERMZY', 'PERMZ',
    'ACTNUM',
    'INCLUDE',
    
]

KeyWordsDatatypes=[#Corrsponding data types
    int,
    int,
    int,int,int,int,
    float,float,
    float,
    float,float,float,
    float,float,float,
    float,float,float,
    int
]

def parseDataArray(DataArray):
        """Parse special dataArray format in GRDECL 
        example:
            5*3.0=[3.0 3.0 3.0 3.0 3.0]
            1.0 2*3.0 5.0=[1.0 3.0 3.0 5.0]
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """

        data=[]
        error_count=0
        for value in DataArray:
            if(is_number(value)==2):
                num,val=value.split('*')
                for i in range(int(num)): data.append(val)
            elif(is_number(value)==1):
                data.append(value)
            else:
                error_count+=1
        
        if(error_count>0):
            print(DataArray)
        
        assert error_count==0, '[Error] Can not find any numeric value!'
        
        return data
    
def is_number(s):
    #Determine a string is a number or not
    #Used in [read_GRDECL] [getBlkdata]
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    
    try: #Special format N*val= [val val val ....]
        num, val = s.split('*')
        return 2
    except ValueError:
        pass
 
    return False

## Auxilary functions
def cell_id(i,j,k,nx,ny):
    """
    Get the cell Id given i,j,k indexes. 
        * ---  *  ---  *  --- *
        | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
        * ---  *  ---  *  --- *
        | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
        * ---  *  ---  *  --- *
        | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
        * ---  *  ---  * ---  *
    """
    cell = (nx*j+i)+k*nx*ny

    return cell

def cell_ijk(cell_id,nx,ny):
    """
    Get the cell indexes i,j,k given the cell id
        * ---  *  ---  *  --- *
        | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
        * ---  *  ---  *  --- *
        | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
        * ---  *  ---  *  --- *
        | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
        * ---  *  ---  * ---  *
    """
    k=math.ceil(cell_id/(nx*ny))-1
    j=math.ceil((cell_id-(nx*ny)*k)/nx)-1
    i=math.ceil(cell_id-(nx*ny*k)-nx*j)
    return i,j,k

#Interpolate z on pillars
def interpolate_z_pillar(z,p):
    """
    Obtain the eight coords for a cell
        X,Y coords has to be interpolated from Z
    xy1=xy0+k*z
    Pillar=np.array([[x0 y0 z0],[x1 y1 z1]])
    """
    x = ((p[0,0]-p[1,0])/(p[0,2]-p[1,2]))*(z-p[0,2])+p[0,0]
    y = ((p[0,1]-p[1,1])/(p[0,2]-p[1,2]))*(z-p[0,2])+p[0,1]

    xyz = np.array([x,y,z])
    return xyz

#3D Rotation funtions
def rotation(points,azimuth,dip,plunge):
    assert points.ndim == 2

    azi_rad = np.radians(azimuth)
    dip_rad = np.radians(dip)
    plg_rad = np.radians(plunge)

    ry = np.array([
        [np.cos(plg_rad),0,-np.sin(plg_rad)],
        [0,1,0],
        [np.sin(plg_rad),0,np.cos(plg_rad)],
    ])
    rx = np.array([
        [1,0,0],
        [0,np.cos(dip_rad),np.sin(dip_rad)],
        [0,-np.sin(dip_rad),np.cos(dip_rad)]
    ])
    rz = np.array([
        [np.cos(azi_rad),-np.sin(azi_rad),0],
        [np.sin(azi_rad),np.cos(azi_rad),0],
        [0,0,1]
    ])

    rot = np.matmul(np.matmul(ry,rx),rz)

    rot_points = np.matmul(points,rot)

    return rot_points

def RemoveCommentLines(data,commenter='--'):
    #Remove comment and empty lines
    data_lines=data.strip().split('\n')
    newdata=[]
    for line in data_lines:
        if line.startswith(commenter) or not line.strip():
            # skip comments and blank lines
            continue   
        newdata.append(line)
    return '\n'.join(newdata)


def scanKeyword(data):
    #scan and find the keyword
    #e.g. ['INIT','DX','2500*30.0'] -> ['DX','2500*30.0']
    for key in SupportKeyWords:
        if (key in data) and (data.find(key)!=0):
            return data[data.find(key):-1]
    return data
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
        self.tops = kwargs.pop('tops',None)  #Pending to define
        self.origin = kwargs.pop('origin',None)
        self.azimuth = kwargs.pop('azimuth',0) #counterclockwise
        self.dip = kwargs.pop('dip',0) #clockwise
        self.plunge = kwargs.pop('plunge',0) #clockwise

        #Corner Point Grid
        self.coord = kwargs.pop('coord',None)
        self.zcorn = kwargs.pop('zcorn',None)

        #Petrophysical properties
        self.spatial_data = kwargs.pop('spatial_data',None)
        
        self.skiped_keywords = 0

#####################################################
############## Properties ###########################

    #---Grid Type---
    @property
    def grid_type(self):
        return self._grid_type

    @grid_type.setter
    def grid_type(self,value):
        if value is not None:
            assert isinstance(value,str), f'{type(value)} not accepted. Name must be str'
            assert value in ['cartesian','corner_point']
        self._grid_type = value

    #---NX---
    @property
    def nx(self):
        return self._nx

    @nx.setter 
    def nx(self,value):
        if value is not None:
            assert isinstance(value,(int,np.int64)) and value > 0, f'{type(value)} not accepted. Must be int'
        self._nx = value 

    #---NY---
    @property
    def ny(self):
        return self._ny

    @ny.setter 
    def ny(self,value):
        if value is not None:
            assert isinstance(value,(int,np.int64)) and value > 0 
        self._ny = value 

    #---NZ---
    @property
    def nz(self):
        return self._nz

    @nz.setter 
    def nz(self,value):
        if value is not None:
            assert isinstance(value,(int,np.int64)) and value > 0 
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
        if value is not None:
            assert isinstance(value,(int,float,list,np.ndarray)), 'Value entered is not the right type Tyepes allowed: (int,float,list,np.ndarray,None)'
            
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
        else:
            self._dx = value

    #---DY---
    @property
    def dy(self):
        assert len(self._dx) == self.n
        return self._dy 

    @dy.setter 
    def dy(self,value):
        if value is not None:
            assert isinstance(value,(int,float,list,np.ndarray)), 'Value entered is not the right type Tyepes allowed: (int,float,list,np.ndarray,None)'
            
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
        else:
            self._dy = value
            
    #---DZ---
    @property
    def dz(self):
        assert len(self._dx) == self.n
        return self._dz 

    @dz.setter 
    def dz(self,value):
        if value is not None:
            assert isinstance(value,(int,float,list,np.ndarray)), 'Value entered is not the right type Tyepes allowed: (int,float,list,np.ndarray,None)'
            
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
        else:
            self._dx = value

    @property
    def origin(self):
        return self._origin

    @origin.setter
    def origin(self,value):
        if value is not None:
            assert isinstance(value,Point), 'Origin point must be Point or None types'
            if value is None:
                assert self._grid_type != 'cartesian', 'If cartesian grid set, origin must be set'
            else:
                assert value.has_z, 'Point must have x,y,z coordinates'
        self._origin = value

    @property
    def azimuth(self):
        return self._azimuth

    @azimuth.setter
    def azimuth(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), 'Must be a number'

            if isinstance(value,np.ndarray):
                _az = value.flatten(order='F')
                assert len(_az) == 1, f'list must length of 1'
            
            assert value >= -360 and value <= 360, 'Azimuth angle must be between 0 and 360'
        self._azimuth = value

    @property
    def dip(self):
        return self._dip

    @dip.setter
    def dip(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), 'Must be a number'

            if isinstance(value,np.ndarray):
                _az = value.flatten(order='F')
                assert len(_az) == 1, f'list must length of 1'
            
            assert value >= -360 and value <= 360, 'dip angle must be between 0 and 360'
        self._dip = value

    @property
    def plunge(self):
        return self._plunge

    @plunge.setter
    def plunge(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), 'Must be a number'

            if isinstance(value,np.ndarray):
                _az = value.flatten(order='F')
                assert len(_az) == 1, f'list must length of 1'
            
            assert value >= -360 and value <= 360, 'plunge angle must be between 0 and 360'
        self._plunge = value

    @property
    def cartesian_vertices_coord(self):
        #Vertices coordinates starting at 0,0,0
        x_vert_cord = np.concatenate((np.zeros(1),self.dx.reshape((self.nx,self.ny,self.nz),order='f')[:,0,0]),axis=0).cumsum()
        y_vert_cord = np.concatenate((np.zeros(1),self.dy.reshape((self.nx,self.ny,self.nz),order='f')[0,:,0]),axis=0).cumsum()
        z_vert_cord = -np.concatenate((np.zeros(1),self.dz.reshape((self.nx,self.ny,self.nz),order='f')[0,0,:]),axis=0).cumsum()

        points = np.zeros(((self.nx+1)*(self.ny+1)*(self.nz+1),3))
        for k in range(self.nz+1):
            for j in range(self.ny+1):
                for i in range(self.nx+1):
                    l = cell_id(i,j,k,self.nx+1,self.ny+1)
                    points[l,0] = x_vert_cord[i]
                    points[l,1] = y_vert_cord[j]
                    points[l,2] = z_vert_cord[k]

        #Get rotated points with respect 0,0,0
        rot_points = rotation(points,self.azimuth,self.dip,self.plunge)
        
        #Adjust the coordinates according with Origin Point
        origin = np.array([self.origin.x,self.origin.y,self.origin.z])

        self._vertices_coord = rot_points + origin
        return self._vertices_coord

    @property
    def cartesian_center_point_coord(self):
        #Vertices coordinates starting at 0,0,0
        x_vert_cord = np.concatenate((np.zeros(1),self.dx.reshape((self.nx,self.ny,self.nz),order='f')[:,0,0]),axis=0).cumsum()
        y_vert_cord = np.concatenate((np.zeros(1),self.dy.reshape((self.nx,self.ny,self.nz),order='f')[0,:,0]),axis=0).cumsum()
        z_vert_cord = -np.concatenate((np.zeros(1),self.dz.reshape((self.nx,self.ny,self.nz),order='f')[0,0,:]),axis=0).cumsum()

        center = np.zeros(((self.nx)*(self.ny)*(self.nz),3))
        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    l = cell_id(i,j,k,self.nx,self.ny)
                    center[l,0] = np.mean((x_vert_cord[i], x_vert_cord[i+1]))
                    center[l,1] = np.mean((y_vert_cord[j], y_vert_cord[j+1]))
                    center[l,2] = np.mean((z_vert_cord[k], z_vert_cord[k+1]))

        #Get rotated points with respect 0,0,0
        rot_points = rotation(center,self.azimuth,self.dip,self.plunge)
        
        #Adjust the coordinates according with Origin Point
        origin = np.array([self.origin.x,self.origin.y,self.origin.z])

        self._center_coord = rot_points + origin
        return self._center_coord


    # Corner Point Grid
    #---COORD---
    @property
    def coord(self):
        return self._coord 

    @coord.setter 
    def coord(self,value):
        if value is not None:
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
        else:
            self._coord = value

    #---ZCORN---
    @property
    def zcorn(self):
        return self._zcorn 

    @zcorn.setter 
    def zcorn(self,value):
        if value is not None:
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
        else:
            self._zcorn = value
    @property 
    def spatial_data(self):
        return self._spatial_data 

    @spatial_data.setter 
    def spatial_data(self,value):
        if value is not None:
            assert isinstance(value,dict)

            if isinstance(value,dict):
                for i in value:
                    #assert i in petrophysical_properties, f"Keyword {i} not in supported properties {petrophysical_properties} "
                    assert isinstance(value[i],(int,float,list,np.ndarray))     

                    if isinstance(value[i],(int,float)):
                        self._spatial_data[i] = np.full(self.n,value[i])
                    elif isinstance(value[i],list):
                        _prop = np.array(value[i]).flatten(order='F')
                        assert len(_prop) == self.n, f'{i} list must be of length {self.n}'
                        assert np.issubdtype(_prop.dtype, np.number), f'{i} List must contain only numbers'
                        assert all(_prop>=0), f'{i} must be greater than 0'
                        self._spatial_data[i] = _prop
                    elif isinstance(value[i],np.ndarray):
                        _prop = value[i].flatten(order='F')
                        assert len(_prop) == self.n, f'{i} list must be of length {self.n}'
                        assert np.issubdtype(_prop.dtype, np.number), f'{i} List must contain only numbers'
                        assert all(_prop>=0), f'{i} must be greater than 0'
                        self._spatial_data[i] = _prop
        else:
            self._spatial_data = {}


#####################################################
############## Methods ###########################

    def add_spatial_data(self,key,array):
        array = np.atleast_1d(array).flatten(order='F')
        assert self.n == array.shape[0]
        try:
            self._spatial_data.update({key:array})
        except Exception as e:
            print(e)
            raise Exception
        else:
            print(f'Added {key}')

    def read_IncludeFile(self,filename_include,NumData):
        """Read Include data file
        this data file just a series of values
        e.g. 0.2 0.3 12.23 ....
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Aug. 2018
        """

        f=open(filename_include)
        contents=f.read()
        block_dataset=contents.strip().split() #Sepeart input file by slash /
        block_dataset=np.array(block_dataset,dtype=float)
        if(len(block_dataset)!=NumData):
            print('Data size %s is not equal to defined block dimension (NX*NY*NZ) %s'%(len(block_dataset),NumData))
        return block_dataset
    
    def LoadVar(self,Keyword,DataArray,DataSize):
        """Load varables into class
        example:
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        if(Keyword in SupportKeyWords):#KeyWords Check
            assert len(DataArray)==DataSize,'\n     [Error-%s] Incompatible data size! %d-%d' %(Keyword,len(DataArray),DataSize)
            KeywordID=SupportKeyWords.index(Keyword)
            print('     [%s] '%(Keyword),end='')
            self._spatial_data[Keyword]=np.array(DataArray,dtype=KeyWordsDatatypes[KeywordID])
        else:
            print('     [Warnning] Unsupport keywords[%s]' % (Keyword))
            self.skiped_keywords+=1
    
    def read_GRDECL(self,file):
        """Read input file(GRDECL) of Reservoir Simulator- Petrel (Eclipse)  
        file format:http://petrofaq.org/wiki/Eclipse_Input_Data
        
        Arguments
        ---------
        NX, NY, NZ -- Grid dimension.
        blockData_raw -- [0] Keywords [1] values
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2017
        """
        debug=0

        print('[Input] Reading ECLIPSE/PETREL file \"%s\" ....'%(file))

        #Read whole file into list
        f=open(file)
        contents=f.read()
        contents=RemoveCommentLines(contents,commenter='--')
        contents_in_block=contents.strip().split('/') #Sepeart input file by slash /
        contents_in_block = [x for x in contents_in_block if x]#Remove empty block at the end
        NumKeywords=len(contents_in_block)
        print(f'Num Keywords {NumKeywords}')
        GoodFlag=0
        for i,block in enumerate(contents_in_block):#Keyword, Block-wise
            #Clean the data where no spliter \ provided
            block=scanKeyword(block)

            blockData_raw=block.strip().split()
            Keyword=''
            DataArray=[]
            if(len(blockData_raw)>1):
                if(blockData_raw[0]=='ECHO'): #This keyword may next to real keyword
                    Keyword,DataArray=blockData_raw[1],blockData_raw[2:]                    
                else:
                    Keyword,DataArray=blockData_raw[0],blockData_raw[1:]

            #Read Grid Dimension [SPECGRID] or [DIMENS] 
            print(Keyword)
            if(Keyword=='DIMENS'):
                DataArray=np.array(DataArray[:3],dtype=int)
                self.grid_type='cartesian'
                self.nx,self.ny,self.nz=DataArray[0],DataArray[1],DataArray[2]
                print("     Grid Type=%s Grid" %(self.grid_type))
                print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.nx,self.ny,self.nz))
                print("     NumOfGrids=%s"%(self.n))
                print('     NumOfKeywords=%s'%(NumKeywords))
                print("     Reading Keyword %d [%s] " %(i+1,Keyword),end='')
                GoodFlag=1
                continue
            elif(Keyword=='SPECGRID'):
                DataArray=np.array(DataArray[:3],dtype=int)
                self.grid_type='corner_point'
                self.nx,self.ny,self.nz=DataArray[0],DataArray[1],DataArray[2]
                print("     Grid Type=%s" %(self.grid_type))
                print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.nx,self.ny,self.nz))
                print("     NumOfGrids=%s"%(self.n))
                print('     NumOfKeywords=%s'%(NumKeywords))
                print("     Reading Keywords [%s] " %(Keyword),end='')
                GoodFlag=1
                continue
            
            if(self.grid_type is None):#Skip unnecessary keywords
                continue

            if(Keyword in SupportKeyWords): #We need parse the special format in 
                if Keyword == 'INCLUDE':
                #if(len(DataArray)==1 and '.' in DataArray[0]):
                    folder_name=os.path.dirname(file)
                    self.read_GRDECL(os.path.join(folder_name,DataArray[0].replace("'","")))
                    continue
                    #DataArray=self.read_IncludeFile(os.path.join(folder_name,DataArray[0]),self.n)
                print(f'------{Keyword}------')

                DataArray=parseDataArray(DataArray)
            

                #Read Grid spatial information, x,y,z ordering
                #Corner point cell
                if(Keyword=='COORD'):# Pillar coords
                    assert len(DataArray)==6*(self.nx+1)*(self.ny+1),'[Error] Incompatible COORD data size!'
                    self.coord=np.array(DataArray,dtype=float)       
                elif(Keyword=='ZCORN'):# Depth coords
                    assert len(DataArray)==8*self.n, '[Error] Incompatible ZCORN data size!'
                    self.zcorn=np.array(DataArray,dtype=float)
                
                #Cartesian cell
                elif(Keyword=='DX'):# Grid size in X dir
                    assert len(DataArray)==self.n, '[Error] Incompatible DX data size!'
                    self.dx=np.array(DataArray,dtype=float)
                elif(Keyword=='DY'):# Grid size in Y dir
                    assert len(DataArray)==self.n, '[Error] Incompatible DY data size!'
                    self.dy=np.array(DataArray,dtype=float)
                elif(Keyword=='DZ'):# Grid size in Z dir
                    assert len(DataArray)==self.n, '[Error] Incompatible DZ data size!'
                    self.dz=np.array(DataArray,dtype=float)
                elif(Keyword=='TOPS'):# TOP position
                    assert len(DataArray)==self.n, '[Error] Incompatible TOPS data size!'
                    self.tops=np.array(DataArray,dtype=float)

                #Read Grid Properties information
                else:
                    self.LoadVar(Keyword,DataArray,DataSize=self.n)

        f.close()
        #assert GoodFlag==1,'Can not find grid dimension info, [SPECGRID] or [DIMENS]!'
        print('.....Done!')


        #Genetrate TOPS for cartesian grid if TOPS if not given
        if(self.grid_type=='Cartesian' and len(self.tops)==0):
            self.tops=np.zeros(self.n)
            for k in range(self.nz-1):
                for j in range(self.ny):
                    for i in range(self.nx):
                        ijk=cell_id(i,j,k,self.nx,self.ny)
                        ijk_next=cell_id(i,j,k+1,self.nx,self.ny)
                        self.tops[ijk_next] = self.tops[ijk] + self.dz[ijk]


    def to_ecl(self, filename=None):
        string = "-- Data Exported from reservoirpy Python Package\n"
        if self.grid_type == 'cartesian':
            
            string += 'TOPS\n'
            string += ' ' + ' '.join([str(v) + '\n' if (i+1)%5==0 else str(v) for i,v in enumerate(self.tops)]) + '/\n'
            
            string += 'DX\n'
            string += ' ' + ' '.join([str(v) + '\n' if (i+1)%5==0 else str(v) for i,v in enumerate(self.dx)]) + '/\n'

            string += 'DY\n'
            string += ' ' + ' '.join([str(v) + '\n' if (i+1)%5==0 else str(v) for i,v in enumerate(self.dy)]) + '/\n'
        
            string += 'DZ\n'
            string += ' ' + ' '.join([str(v) + '\n' if (i+1)%5==0 else str(v) for i,v in enumerate(self.dz)]) + '/\n'
        
        # elif corner_point
        else:
            
            string += 'SPECGRID\n'
            string += f' {self.nx} {self.ny} {self.nz} 1 F\n'                               

            print('COORD')
            string += 'COORD\n'
            #string += ' ' + pd.DataFrame(self.coord.reshape(-1,6)).to_string(index=False,header=False) + '\n'
            string += ' ' + ' '.join([str(v) + '\n' if (i+1)%10==0 else str(v) for i,v in enumerate(self.coord)]) + '/\n'
            print('ZCOORN')
            string += 'ZCORN\n'
            #string += ' ' + pd.DataFrame(self.zcorn.reshape(-1,8)).to_string(index=False,header=False) + '\n'
            string += ' ' + ' '.join([str(v) + '\n' if (i+1)%10==0 else str(v) for i,v in enumerate(self.zcorn)]) + '/\n'
        if bool(self.spatial_data):
            for key in self.spatial_data.keys():
                print(key)
                string += key + '\n'
                string += ' ' + ' '.join([str(v) + '\n' if (i+1)%10==0 else str(v) for i,v in enumerate(self.spatial_data[key])])  + '/\n'
 
        if filename is not None:
            try:
                with open(filename,'w') as text_file:
                    text_file.write(string)
            except Exception as e:
                print(e)
                pass
        return string
            
    
    def get_cell_id(self,i,j,k):
        """
        Get the cell Id given i,j,k indexes. 
            * ---  *  ---  *  --- *
            | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
            * ---  *  ---  *  --- *
            | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
            * ---  *  ---  *  --- *
            | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
            * ---  *  ---  * ---  *
        """
        c_id = cell_id(i,j,k,self.nx,self.ny)
        return c_id

    def get_cell_ijk(self,cell_id):
        """
        Get the cell indexes i,j,k given the cell id
            * ---  *  ---  *  --- *
            | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
            * ---  *  ---  *  --- *
            | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
            * ---  *  ---  *  --- *
            | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
            * ---  *  ---  * ---  *
        """
        i,j,k=cell_ijk(cell_id,self.nx,self.ny)
        return i,j,k

    def get_pillar(self,pillar_id:int):
        """
        Get the Top and Bottom coordinates of a pillar id
        """
        if self.grid_type == 'corner_point':
            id_top=[6*pillar_id+0,6*pillar_id+1,6*pillar_id+2]
            id_bottom=[6*pillar_id+3,6*pillar_id+4,6*pillar_id+5]
            top_point=np.array([self.coord[i] for i in id_top])
            bottom_point=np.array([self.coord[i] for i in id_bottom])
        else:
            raise ValueError('Pillar are only set in a Corner Point Grid')
        return np.array([top_point,bottom_point])

    def get_cell_pillars(self,i,j):
        """Obtain the four pillars (p0,p1,p2,p3) of a corner point cell
        The index of pillar
        
        3x3x1 system (2D X-Y plane)
        12--- 13  --- 14  ---15
        |      |       |      |  <- Cell 6,7,8
        8 ---  9  --- 10  ---11
        |      |       |      |  <- Cell 3,4,5
        4 ---  5  ---  6  --- 7
        |      |       |      |  <- Cell 0,1,2
        0 ---  1 ---   2 ---  3
        
        The pillars index for a grid follows below ordering (XY Plane)
        p2   p3
        *------*
        |      |
        |      |
        *------*
        p0   p1

        """
        if self.grid_type == 'corner_point':
            p0 = cell_id(i,j,0,self.nx+1,self.ny+1)
            p1 = cell_id(i+1,j,0,self.nx+1,self.ny+1)
            p2 = cell_id(i,j+1,0,self.nx+1,self.ny+1)
            p3 = cell_id(i+1,j+1,0,self.nx+1,self.ny+1)

            pls = [self.get_pillar(p0),self.get_pillar(p1),self.get_pillar(p2),self.get_pillar(p3)]
        else:
            raise ValueError('Pillar are only set in a Corner Point Grid')
        return np.array(pls)

    def get_vertices_id(self,i,j,k,order='GRD'):
        """
        Cartesian Grid

        Get the cell Id given i,j,k indexes. 
            13 --- 14  --- 15 --- 16
            | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
            9 ---  10  --- 11 --- 12
            | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
            5 ---  6  ---  7  --- 8
            | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
            1 ---  2  ---  3 ---  4

        Corner Point Grid

        3x3x1 system (2D X-Y plane)
        30---31,32---33,34---35
        |      |       |      |  <- Cell 6,7,8
        24---25,26---27,28---29
        18---19,20---21,22---23
        |      |       |      |  <- Cell 3,4,5
        12---13,14---15,16---17
        6 --- 7,8 --- 9,10---11
        |      |       |      |  <- Cell 0,1,2
        0 --- 1,2 --- 3,4 --- 5
        Node order convention for a 3D cell
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1
        """ 
        if self.grid_type == 'cartesian':
            nx,ny=self.nx+1,self.ny+1
            p0=cell_id(i,j,k,nx,ny)
            p1=cell_id(i+1,j,k,nx,ny)
            p2=cell_id(i,j+1,k,nx,ny)
            p3=cell_id(i+1,j+1,k,nx,ny)

            p4=cell_id(i,j,k+1,nx,ny)
            p5=cell_id(i+1,j,k+1,nx,ny)
            p6=cell_id(i,j+1,k+1,nx,ny)
            p7=cell_id(i+1,j+1,k+1,nx,ny)

            if order == 'GRD':
                points = [p0,p1,p2,p3,p4,p5,p6,p7]
            elif order == 'VTK':
                points = [p4,p5,p7,p6,p0,p1,p3,p2]

            return np.array(points)

        if self.grid_type == 'corner_point':
            nx,ny=2*self.nx,2*self.ny
            p0=cell_id(2*i,2*j,2*k,nx,ny)
            p1=cell_id(2*i+1,2*j,2*k,nx,ny)
            p2=cell_id(2*i,2*j+1,2*k,nx,ny)
            p3=cell_id(2*i+1,2*j+1,2*k,nx,ny)

            p4=cell_id(2*i,2*j,2*k+1,nx,ny)
            p5=cell_id(2*i+1,2*j,2*k+1,nx,ny)
            p6=cell_id(2*i,2*j+1,2*k+1,nx,ny)
            p7=cell_id(2*i+1,2*j+1,2*k+1,nx,ny)

            if order == 'GRD':
                points = [p0,p1,p2,p3,p4,p5,p6,p7]
            elif order == 'VTK':
                points = [p4,p5,p7,p6,p0,p1,p3,p2]

            return np.array(points)


    def get_vertices_z(self,i,j,k):
        """
        Node order convention for a 3D cell
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1
        """
        # Get the z coord for a cell
        if self.grid_type == 'corner_point':
            p = self.get_vertices_id(i,j,k)
            z = [self.zcorn[i] for i in p]
            return np.array(z)
        elif self.grid_type == 'cartesian':
            p= self.get_vertices_id(i,j,k)
            z = [self.cartesian_vertices_coord[i,2] for i in p]
            return np.array(z)

        # Pending for cartessian grid

    def get_vertices_coords(self,i,j,k,order='GRD'):
        if self.grid_type == 'corner_point':
            coords=[]
            pillars = self.get_cell_pillars(i,j)
            cell_z =self.get_vertices_z(i,j,k)

            for i in range(8):
                p_id = i%4
                coords.append(interpolate_z_pillar(cell_z[i],pillars[p_id]))
            
            if order == 'GRD':
                v_coord = np.array(coords)
            elif order == 'VTK':
                v_coord = np.array(coords)[[4,5,7,6,0,1,3,2],:]

            return v_coord

        elif self.grid_type == 'cartesian':
            p= self.get_vertices_id(i,j,k)
            coords = [[self.cartesian_vertices_coord[i,0],self.cartesian_vertices_coord[i,1],self.cartesian_vertices_coord[i,2]] for i in p]
            if order == 'GRD':
                v_coord = np.array(coords)
            elif order == 'VTK':
                v_coord = np.array(coords)[[4,5,7,6,0,1,3,2],:]

            return v_coord


    def get_vertices_face_z(self,i,j,k, face=None):
        """
         Get the Z coords for a cell
        
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1   
        Follow getCornerPointCellIdx convention:
        X-, [0,2,4,6]
        X+, [1,3,5,7]
        Y-, [0,1,4,5]
        Y+, [2,3,6,7]
        Z+, [0,1,2,3]
        Z-, [4,5,6,7]
        """
        assert face is not None, 'A face must be choosen'
        points_id = self.get_vertices_id(i,j,k)
        if(face=="X-"): face_id=[points_id[0],points_id[2],points_id[4],points_id[6]]
        if(face=="X+"): face_id=[points_id[1],points_id[3],points_id[5],points_id[7]]
        if(face=="Y-"): face_id=[points_id[0],points_id[1],points_id[4],points_id[5]]
        if(face=="Y+"): face_id=[points_id[2],points_id[3],points_id[6],points_id[7]]
        if(face=="Z-"): face_id=[points_id[4],points_id[5],points_id[6],points_id[7]]
        if(face=="Z+"): face_id=[points_id[0],points_id[1],points_id[2],points_id[3]]

        if self.grid_type == 'cartesian':
            z_face = [self.cartesian_vertices_coord[i,2] for i in face_id]
        elif self.grid_type == 'corner_point':
            z_face = [self.zcorn[i] for i in face_id]

        return np.array(z_face)

    def get_vertices_face_coords(self,i,j,k, face=None):
        """
         Get the Z coords for a cell
        
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1   
        Follow getCornerPointCellIdx convention:
        X-, [0,2,4,6]
        X+, [1,3,5,7]
        Y-, [0,1,4,5]
        Y+, [2,3,6,7]
        Z+, [0,1,2,3]
        Z-, [4,5,6,7]
        """
        assert face is not None, 'A face must be choosen'
        points_id = self.get_vertices_id(i,j,k)
        if (face=="X-"): 
            face_id=[points_id[0],points_id[2],points_id[4],points_id[6]]
            ind = [0,2,4,6]
        elif (face=="X+"): 
            face_id=[points_id[1],points_id[3],points_id[5],points_id[7]]
            ind = [1,3,5,7]
        elif (face=="Y-"): 
            face_id=[points_id[0],points_id[1],points_id[4],points_id[5]]
            ind = [0,1,4,5]
        elif (face=="Y+"): 
            face_id=[points_id[2],points_id[3],points_id[6],points_id[7]]
            ind = [2,3,6,7]
        elif (face=="Z-"): 
            face_id=[points_id[4],points_id[5],points_id[6],points_id[7]]
            ind = [4,5,6,7]
        elif (face=="Z+"): 
            face_id=[points_id[0],points_id[1],points_id[2],points_id[3]]
            ind = [0,1,2,3]

        if self.grid_type == 'cartesian':
            z_face = [self.cartesian_vertices_coord[i,:] for i in face_id]
        elif self.grid_type == 'corner_point':
            v_cord = self.get_vertices_coords(i,j,k)
            z_face = v_cord[ind,:]

        return np.array(z_face)


    def get_center_coord(self,i,j,k):
        """
        Get the cell Id given i,j,k indexes. 
            * ---  *  ---  *  --- *
            | 0,2  |  1,2  |  2,2 |  <- Cell id 6,7,8
            * ---  *  ---  *  --- *
            | 0,1  |  1,1  |  2,1 |  <- Cell id 3,4,5
            * ---  *  ---  *  --- *
            | 0,0  |  1,0  |  2,0 |  <- Cell id 0,1,2
            * ---  *  ---  * ---  *
        """
        cid = self.get_cell_id(i,j,k)
        if self.grid_type == 'cartesian':
            center = self.cartesian_center_point_coord[cid,:]
        elif self.grid_type == 'corner_point':
            points = self.get_vertices_coords(i,j,k)
            center = points.mean(axis=0)
        
        return center

    def get_vtk(self):
        """
        Get the pyvista Object
        https://docs.pyvista.org/examples/00-load/create-unstructured-surface.html#sphx-glr-examples-00-load-create-unstructured-surface-py
        """
 
        #Identify the cell data connections
        offset = np.arange(0,9*self.n,step=9)

        points = np.zeros((self.n*8,3))
        #Cells
        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    c_idx = self.get_cell_id(i,j,k)
                    #cells_array[c_idx,:] = self.get_vertices_id(i,j,k, order='VTK')

                    ind_from = 8*c_idx
                    ind_to = 8*(c_idx+1)
                    points[ind_from:ind_to,:] = self.get_vertices_coords(i,j,k, order = 'VTK')

        # Make a vector of shape self.n, make 2D and append to cell array then flatten C order
        cell_array = np.arange(self.n*8).reshape((self.n,8))
        cells = np.append(np.full(self.n,8).reshape((self.n,1)),cell_array,1).flatten()

        # cell type array. Contains the cell type of each cell
        cell_type = np.array([vtk.VTK_HEXAHEDRON]*self.n)

        grid = pv.UnstructuredGrid(offset, cells, cell_type, points)

        if self.spatial_data is not None:
            for i in self.spatial_data.items():
                grid.cell_arrays[i[0]] = i[1]

        return grid







        


                






    
    

    

    




