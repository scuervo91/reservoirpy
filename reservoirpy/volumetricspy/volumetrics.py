import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import pandas as pd
from skimage import measure
from scipy.integrate import simps

def poly_area(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

class surface:
    def __init__(self, **kwargs):
        self.x = kwargs.pop('x',None)
        self.y = kwargs.pop('y',None)
        self.z = kwargs.pop('z',None)

    #Properties

    @property
    def x(self):
        return self._x 

    @x.setter 
    def x(self,value):
        if value is not None:
            assert isinstance(value,np.ndarray)
            assert value.ndim == 2 
        self._x = value

    @property
    def y(self):
        return self._y 

    @y.setter 
    def y(self,value):
        if value is not None:
            assert isinstance(value,np.ndarray)
            assert value.ndim == 2 
        self._y = value

    @property
    def z(self):
        return self._z 

    @z.setter 
    def z(self,value):
        if value is not None:
            assert isinstance(value,np.ndarray)
            assert value.ndim == 2 
        self._z = value

    def contour(self,ax=None,**kwargs):

        #Create the Axex
        cax= ax or plt.gca()

        cax.contour(self.x,self.y,self.z,**kwargs)

    def contourf(self,ax=None,**kwargs):

        #Create the Axex
        cax= ax or plt.gca()

        cax.contourf(self.x,self.y,self.z,**kwargs)
    
    def structured_surface_vtk(self):

        #Get a Pyvista Object StructedGrid
        grid = pv.StructuredGrid(self.x, self.y, self.z).elevation()

        return grid

    
    def get_contours(self,levels=None,n=10):
        
        #define levels
        if levels is not None:
            assert isinstance(levels,(np.ndarray,list))
            levels = np.atleast_1d(levels)
            assert levels.ndim==1
        else:
            zmin = np.nanmin(self.z)
            zmax = np.nanmax(self.z)

            levels = np.linspace(zmin,zmax,n)

        zz = self.z
        xmax = np.nanmax(self.x)
        ymax = np.nanmax(self.y)
        xmin = np.nanmin(self.x)
        ymin = np.nanmin(self.y)

        #iterate over levels levels
        data = pd.DataFrame()
        i = 0
        for level in levels:
            contours = measure.find_contours(zz,level)

            if contours == []:
                continue
            else:
                for contour in contours:
                    level_df = pd.DataFrame(contour, columns=['y','x'])
                    level_df['level'] = level
                    level_df['n'] = i
                    data = data.append(level_df,ignore_index=True)
                    i += 1

        if not data.empty:
            #re scale
            data['x'] = (data['x']/zz.shape[1]) * (xmax - xmin) + xmin
            data['y'] = (data['y']/zz.shape[0]) * (ymax - ymin) + ymin

        return data

    def get_contours_area(self,levels=None,n=10, group=True,c=2.4697887e-4):

        #c is the conversion factor from m2 to acre

        #get contours
        contours = self.get_contours(levels=levels,n=n)

        if contours.empty:
            raise ValueError('None contours found')
        #dataframe
        data = pd.DataFrame()

        for level in contours['level'].unique():
            level_df = contours.loc[contours['level']==level,:]

            for n in level_df['n'].unique():
                poly_df = level_df.loc[level_df['n']==n,:]

                area = poly_area(poly_df['x'].values, poly_df['y'].values) * c

                area_df = pd.DataFrame({'level':[level],'n':[n],'area':area})

                data = data.append(area_df,ignore_index=True)

        if group:
            data_g = data[['level','area']].groupby('level').sum()
            return data_g
        else:
            return data

    def get_volume(self,levels=None, n=10,c=2.4697887e-4):
        
        
        area = self.get_contours_area(levels=levels,n=n,c=c,group=True)

        #Integrate
        rv=simps(np.abs(area.index),area['area'])

        return rv, area

class surface_group:
    def __init__(self,**kwargs):
       
        self.surfaces = kwargs.pop('surfaces',None) 

    @property
    def surfaces(self):
        return self._surfaces

    @surfaces.setter 
    def surfaces(self,value):
        if value is not None:
            assert isinstance(value,dict)
            assert all(isinstance(value[i],surface) for i in value)
            self._surfaces = value
        else:
            self._surfaces = {}

    def add_surface(self,surf):
        assert isinstance(surf,dict)
        assert all(isinstance(surf[i],surface) for i in surf)

        _surface_dict = self.surfaces.copy()

        _surface_dict.update(surf)
        self._surfaces = _surface_dict

    def get_volume(self, top_surface=None, bottom_surface=None, levels=None, n=10,c=2.4697887e-4):

        assert all([top_surface is not None,bottom_surface is not None])

        top_area = self.surfaces[top_surface].get_contours_area(levels=levels,n=n,c=c,group=True)
        bottom_area = self.surfaces[bottom_surface].get_contours_area(levels=levels,n=n,c=c,group=True)

        #Merge two contours ara for top and bottom indexed by depth
        area=top_area.merge(bottom_area,how='outer',left_index=True,right_index=True,suffixes=['_top','_bottom']).fillna(0)
        area['dif']= area['area_top'] - area['area_bottom']

        #Integrate
        rv=simps(np.abs(area.index),area['dif'])

        return rv, area








            
            





    

