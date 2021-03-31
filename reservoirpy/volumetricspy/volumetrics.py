import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import pandas as pd
from skimage import measure
from scipy.integrate import simps
from scipy.interpolate import griddata
import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon
from zmapio import ZMAPGrid

def poly_area(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

class Surface:
    def __init__(self, **kwargs):
        self.x = kwargs.pop('x',None)
        self.y = kwargs.pop('y',None)
        self.z = kwargs.pop('z',None)
        self.crs = kwargs.pop('crs',4326)

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

    @property
    def crs(self):
        return self._crs

    @crs.setter
    def crs(self,value):
        assert isinstance(value,(int,str,type(None))), f"{type(value)} not accepted. Name must be str. Example 'EPSG:3117'"
        
        if isinstance(value,int):
            value = f'EPSG:{value}'
        elif isinstance(value,str):
            assert value.startswith('EPSG:'), 'if crs is string must starts with EPSG:. If integer must be the Coordinate system reference number EPSG http://epsg.io/'
        self._crs = value

    def contour(self,ax=None,**kwargs):

        #Create the Axex
        cax= ax or plt.gca()

        return cax.contour(self.x,self.y,self.z,**kwargs)

    def contourf(self,ax=None,**kwargs):

        #Create the Axex
        cax= ax or plt.gca()

        return cax.contourf(self.x,self.y,self.z,**kwargs)
        
    
    def structured_surface_vtk(self):

        #Get a Pyvista Object StructedGrid
        grid = pv.StructuredGrid(self.x, self.y, self.z).elevation()

        return grid

    def get_contours_bound(self,levels=None,zmin=None,zmax=None,n=10):
        #define levels
        if levels is not None:
            assert isinstance(levels,(np.ndarray,list))
            levels = np.atleast_1d(levels)
            assert levels.ndim==1
        else:
            zmin = zmin if zmin is not None else np.nanmin(self.z)
            zmax = zmax if zmax is not None else np.nanmax(self.z)

            levels = np.linspace(zmin,zmax,n)

        xmax = np.nanmax(self.x)
        ymax = np.nanmax(self.y)
        xmin = np.nanmin(self.x)
        ymin = np.nanmin(self.y)

        #iterate over levels levels
        contours = self.structured_surface_vtk().contour(isosurfaces=levels.tolist())

        contours.points[:,2] = contours['Elevation']

        df = pd.DataFrame(contours.points, columns=['x','y','z'])

        #Organize the points according their angle with respect the centroid. This is done with the 
        #porpuse of plot the bounds continously.
        list_df_sorted = []

        for i in df['z'].unique():
            df_z = df.loc[df['z']==i,['x','y','z']]
            centroid = df_z[['x','y']].mean(axis=0).values
            df_z[['delta_x','delta_y']] = df_z[['x','y']] - centroid
            df_z['angle'] = np.arctan2(df_z['delta_y'],df_z['delta_x'])
            df_z.sort_values(by='angle', inplace=True)


            list_df_sorted.append(df_z)


        return pd.concat(list_df_sorted, axis=0)


    def get_contours_area_bounds(self,levels=None,n=10,zmin=None,zmax=None,c=2.4697887e-4):

        contours = self.get_contours_bound(levels=levels,zmin=zmin,zmax=zmax,n=n)


        area_dict= {}
        for i in contours['z'].unique():
            poly = contours.loc[contours['z']==i,['x','y']]
            area = poly_area(poly['x'],poly['y'])
            area_dict.update({i:area*c})

        return pd.DataFrame.from_dict(area_dict, orient='index', columns=['area'])

    def get_contours_area_mesh(self,levels=None,n=10,zmin=None,zmax=None,c=2.4697887e-4):


        zmin = zmin if zmin is not None else np.nanmin(self.z)
        zmax = zmax if zmax is not None else np.nanmax(self.z)

        if levels is not None:
            assert isinstance(levels,(np.ndarray,list))
            levels = np.atleast_1d(levels)
            assert levels.ndim==1
        else:
            levels = np.linspace(zmin,zmax,n)

        dif_x = np.diff(self.x,axis=1).mean(axis=0)
        dif_y = np.diff(self.y,axis=0).mean(axis=1)
        dxx, dyy = np.meshgrid(dif_x,dif_y) 
        
        area_dict = {}
        for i in levels:
            z = self.z.copy()
            z[(z<i)|(z>zmax)|(z<zmin)] = np.nan
            z = z[1:,1:]
            a = dxx * dyy * ~np.isnan(z) *2.4697887e-4
            area_dict.update({i:a.sum()})

        return pd.DataFrame.from_dict(area_dict, orient='index', columns=['area'])

    def get_contours(self,levels=None,zmin=None,zmax=None,n=10):
        
        #define levels
        if levels is not None:
            assert isinstance(levels,(np.ndarray,list))
            levels = np.atleast_1d(levels)
            assert levels.ndim==1
        else:
            zmin = zmin if zmin is not None else np.nanmin(self.z)
            zmax = zmax if zmax is not None else np.nanmax(self.z)

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

    def get_contours_gdf(self,levels=None,zmin=None,zmax=None,n=10, crs="EPSG:4326"):
        
        #define levels
        if levels is not None:
            assert isinstance(levels,(np.ndarray,list))
            levels = np.atleast_1d(levels)
            assert levels.ndim==1
        else:
            zmin = zmin if zmin is not None else np.nanmin(self.z)
            zmax = zmax if zmax is not None else np.nanmax(self.z)

            levels = np.linspace(zmin,zmax,n)

        zz = self.z
        xmax = np.nanmax(self.x)
        ymax = np.nanmax(self.y)
        xmin = np.nanmin(self.x)
        ymin = np.nanmin(self.y)

        #iterate over levels levels
        data = gpd.GeoDataFrame()
        i = 0
        for level in levels:
            poly_list =[]
            contours = measure.find_contours(zz,level)

            if contours == []:
                continue
            else:
                for contour in contours:
                    level_df = pd.DataFrame(contour, columns=['y','x'])

                    #Re scale
                    level_df['x'] = (level_df['x']/zz.shape[1]) * (xmax - xmin) + xmin
                    level_df['y'] = (level_df['y']/zz.shape[0]) * (ymax - ymin) + ymin

                    #List of tuples
                    records = level_df[['x','y']].to_records(index=False)
                    list_records = list(records)

                    if len(list_records)<3:
                        continue
                    else:
                        poly = Polygon(list(records))

                    #Append to list of Polygon
                    poly_list.append(poly)

            # Make Multipolygon
            multi_poly = MultiPolygon(poly_list)

            #Make a Geo dataframe
            level_gdf = gpd.GeoDataFrame({'level':[level],'geometry':[multi_poly]})

            # Append data to general geodataframe
            data = data.append(level_gdf,ignore_index=True)
            i += 1
        
        #Add data crs
        data.crs = self.crs 
        
        #Convert to defined crs
        data = data.to_crs(crs)
        return data

    def get_contours_area(self,levels=None,n=10, group=True,c=2.4697887e-4):

        #c is the conversion factor from m2 to acre

        #get contours
        contours = self.get_contours(levels=levels,n=n)

        if contours.empty:
            print('None contours found')
            return pd.Series(np.zeros(levels.shape[0]), index=levels, name='area')
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
        rv=simps(area['area'],np.abs(area.index))

        return rv, area

    def from_z_map(self,value, factor_z = -1, crs=4326):

        z_file = ZMAPGrid(value)
        z_df = z_file.to_dataframe().dropna()
        z_df['Z'] *= factor_z
        p = z_df.pivot(index='Y',columns='X',values='Z')
        p.sort_index(axis=0, inplace=True)
        p.sort_index(axis=1, inplace=True)
        xx,yy = np.meshgrid(p.columns,p.index)

        self.x = xx
        self.y = yy
        self.z = p.values
        self.crs=crs

    def get_z(self, x, y, method='linear'):

        _x = self.x.flatten()
        _y = self.y.flatten()
        _z = self.z.flatten()

        _xf = _x[~np.isnan(_z)]
        _yf = _y[~np.isnan(_z)]
        _zf = _z[~np.isnan(_z)] 

        return griddata((_xf,_yf),_zf,(x,y), method=method)


class SurfaceGroup:
    def __init__(self,**kwargs):
       
        self.surfaces = kwargs.pop('surfaces',None) 

    @property
    def surfaces(self):
        return self._surfaces

    @surfaces.setter 
    def surfaces(self,value):
        if value is not None:
            assert isinstance(value,dict)
            assert all(isinstance(value[i],Surface) for i in value)
            self._surfaces = value
        else:
            self._surfaces = {}

    def add_surface(self,surf):
        assert isinstance(surf,dict)
        assert all(isinstance(surf[i],Surface) for i in surf)

        _surface_dict = self.surfaces.copy()

        _surface_dict.update(surf)
        self._surfaces = _surface_dict

    def get_volume_bounds(self, 
        top_surface=None, 
        bottom_surface=None, 
        levels=None, 
        zmin=None,
        zmax=None,
        n=20,c=2.4697887e-4,method='mesh'):

        assert all([top_surface is not None,bottom_surface is not None])

        #define levels
        if levels is not None:
            assert isinstance(levels,(np.ndarray,list))
            levels = np.atleast_1d(levels)
            assert levels.ndim==1
        else:
            zmin = zmin if zmin is not None else np.nanmin(self.surfaces[bottom_surface].z)
            zmax = zmax if zmax is not None else np.nanmax(self.surfaces[top_surface].z)

            levels = np.linspace(zmin,zmax,n)

        if method=='mesh':
            top_area = self.surfaces[top_surface].get_contours_area_mesh(levels=levels,n=n,c=c,zmin=zmin, zmax=zmax)
            bottom_area = self.surfaces[bottom_surface].get_contours_area_mesh(levels=levels,n=n,c=c,zmin=zmin, zmax=zmax)

        else:
            top_area = self.surfaces[top_surface].get_contours_area_bounds(levels=levels,n=n,c=c,zmin=zmin, zmax=zmax)
            bottom_area = self.surfaces[bottom_surface].get_contours_area_bounds(levels=levels,n=n,c=c, zmin=zmin, zmax=zmax)

        #Merge two contours ara for top and bottom indexed by depth
        area=top_area.merge(bottom_area,how='outer',left_index=True,right_index=True,suffixes=['_top','_bottom']).fillna(0)
        area['dif_area']= np.abs(area['area_top'] - area['area_bottom'])
        area['height'] = np.diff(area.index, append=0)
        area['vol'] = area['dif_area'].multiply(area['height'])
        rv = area['vol'].iloc[0:-1].sum()
        #area['height'] = area.index-area.index.min()
        #area['tick']=np.diff(area['height'], prepend=0)
        #area['vol'] = area['dif_area'] * area['tick']
        #Integrate
        #rv=simps(area['dif'],area['thick'])
        #rv=area['vol'].sum()

        return rv, area.iloc[0:-1]


    def get_volume(self, 
        top_surface=None, 
        bottom_surface=None, 
        levels=None, 
        zmin=None,
        zmax=None,
        n=20,c=2.4697887e-4):

        assert all([top_surface is not None,bottom_surface is not None])

        #define levels
        if levels is not None:
            assert isinstance(levels,(np.ndarray,list))
            levels = np.atleast_1d(levels)
            assert levels.ndim==1
        else:
            zmin = zmin if zmin is not None else np.nanmin(self.surfaces[bottom_surface].z)
            zmax = zmax if zmax is not None else np.nanmax(self.surfaces[top_surface].z)

            levels = np.linspace(zmin,zmax,n)

        top_area = self.surfaces[top_surface].get_contours_area(levels=levels,n=n,c=c,group=True)
        bottom_area = self.surfaces[bottom_surface].get_contours_area(levels=levels,n=n,c=c,group=True)

        #Merge two contours ara for top and bottom indexed by depth
        area=top_area.merge(bottom_area,how='outer',left_index=True,right_index=True,suffixes=['_top','_bottom']).fillna(0)
        area['dif_area']= np.abs(area['area_top'] - area['area_bottom'])
        area['height'] = area.index-area.index.min()
        area['tick']=np.diff(area['height'], prepend=0)
        area['vol'] = area['dif_area'] * area['tick']
        #Integrate
        #rv=simps(area['dif'],area['thick'])
        rv=area['vol'].sum()

        return rv, area

    def structured_surface_vtk(self, surfaces:list=None):
        
        if surfaces is None:
            _surface_list = []
            for key in self.surfaces:
                _surface_list.append(key)
        else:
            _surface_list = surfaces

        data={}
        for s in _surface_list:
        #Get a Pyvista Object StructedGrid
            data[s] = self.surfaces[s].structured_surface_vtk()

        grid_blocks = pv.MultiBlock(data)

        return grid_blocks







            
            





    

