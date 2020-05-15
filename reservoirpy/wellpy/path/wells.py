import pandas as pd 
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from .mincurve import min_curve_method
from .interpolate import interpolate_deviation, interpolate_position
from scipy.interpolate import interp1d
from scipy.spatial import distance_matrix
from pyproj import Proj, transform
import folium
from ...welllogspy.log import log


class perforations(gpd.GeoDataFrame):

    def __init__(self, *args, **kwargs):                                                                                                                                   
        super(perforations, self).__init__(*args, **kwargs)
    
    def get_tick(self):
        try:
            self['md_tick'] = self['md_bottom'] - self['md_top']
        except:
            pass

        try:
            self['tvd_tick'] = self['tvd_bottom'] - self['tvd_top']
        except:
            pass
        
        return self
    
    def get_mid_point(self):
        try:
            self['md_mid_point'] = (self['md_bottom'] + self['md_top'])*0.5
        except:
            pass

        try:
            self['tvd_mid_point'] = (self['tvd_bottom'] + self['tvd_top'])*0.5
        except:
            pass
    
        try:
            self['tvdss_mid_point'] = (self['tvdss_bottom'] + self['tvdss_top'])*0.5
        except:
            pass
        return self
        
    
    @property
    def _constructor(self):
        return perforations
    
class tops(gpd.GeoDataFrame):

    def __init__(self, *args, **kwargs):                                                                                                                                   
        super(tops, self).__init__(*args, **kwargs)     
    
    @property
    def _constructor(self):
        return tops
  
    
class well:
    def __init__(self, **kwargs):

        self.name = kwargs.pop('name', None)
        self.rte = kwargs.pop('rte', None)
        self.surf_coord = kwargs.pop('surf_coord', None)
        self.crs = kwargs.pop('crs', None)
        self.perforations = kwargs.pop('perforations', None)
        self.tops = kwargs.pop('tops', None)
        self.openlog = kwargs.pop('openlog', None)
        self.masterlog = kwargs.pop('masterlog', None) 
        self.caselog = kwargs.pop('caselog', None)
        self.deviation = kwargs.pop('deviation', None) # pandas md inc azi
        self._survey = None

#####################################################
############## Properties ###########################

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,value):
        assert isinstance(value,(str,type(None))), f'{type(value)} not accepted. Name must be str'
        self._name = value

    @property
    def rte(self):
        return self._rte

    @rte.setter
    def rte(self,value):
        assert isinstance(value,(int,float,type(None))), f'{type(value)} not accepted. Name must be number'
        self._rte = value

    @property
    def surf_coord(self):
        return self._surf_coord

    @surf_coord.setter
    def surf_coord(self,value):
        assert isinstance(value,(Point,type(None))), f'{type(value)} not accepted. Name must be Point'
        self._surf_coord = value

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

    @property
    def perforations(self):
        return self._perforations

    @perforations.setter
    def perforations(self,value):
        assert isinstance(value,(perforations,type(None))), f'{type(value)} not accepted. Name must be reservoirpy.wellpy.path.perforations'
        self._perforations = value

    @property
    def tops(self):
        return self._tops

    @tops.setter
    def tops(self,value):
        assert isinstance(value,(tops,type(None))), f'{type(value)} not accepted. Name must be reservoirpy.wellpy.path.tops'
        self._tops = value    

    @property
    def openlog(self):
        return self._openlog

    @openlog.setter
    def openlog(self,value):
        assert isinstance(value,(log,type(None))), f'{type(value)} not accepted. Name must be reservoirpy.welllogspy.log.log'
        self._openlog = value

    @property
    def masterlog(self):
        return self._masterlog

    @masterlog.setter
    def masterlog(self,value):
        assert isinstance(value,(log,type(None))), f'{type(value)} not accepted. Name must be reservoirpy.welllogspy.log.log'
        self._masterlog = value

    @property
    def caselog(self):
        return self._caselog

    @caselog.setter
    def caselog(self,value):
        assert isinstance(value,(log,type(None))), f'{type(value)} not accepted. Name must be reservoirpy.welllogspy.log.log'
        self._caselog = value

    @property
    def deviation(self):
        return self._deviation

    @deviation.setter
    def deviation(self,value):
        assert isinstance(value,(pd.DataFrame,type(None))), f'{type(value)} not accepted. Name must be pd.DataFrame'
        if isinstance(value,pd.DataFrame):
            assert all(i in ['md','inc','azi'] for i in value.columns), "pd.DataFrame must contain columns ['md','inc','azi']"
        self._deviation = value

    @property
    def survey(self):
        if self._deviation is not None:
            self._survey = min_curve_method(
                self._deviation['md'],
                self._deviation['inc'],
                self._deviation['azi'],
                surface_easting=self._surf_coord.x, 
                surface_northing=self._surf_coord.y, 
                kbe=self._rte,
                crs=self._crs)
        return self._survey

#####################################################
############## methods ###########################

    def sample_deviation(self,step=100):
        if self._deviation is not None:
            _survey = self.survey
            new_dev = interpolate_deviation(_survey.index, 
                                            _survey['inc'], 
                                            _survey['azi'], md_step=step)
        else:
            raise ValueError("No survey has been set")
        return new_dev

    def sample_position(self,step=100):
        if self._deviation is not None:
            _survey = self.survey
            new_pos = interpolate_position(_survey['tvd'], 
                                            _survey['easting'], 
                                            _survey['northing'], 
                                            tvd_step=step)
            new_pos_gpd = gpd.GeoDataFrame(new_pos,geometry=gpd.points_from_xy(new_pos.new_easting,new_pos.new_northing),crs=self._crs)
        else:
            raise ValueError("No survey has been set")
        return new_pos_gpd

    """
            #Set the depth interpolators
            self._tvd_int = interp1d(self.survey.index,self.survey['tvd'])
            self._tvdss_int = interp1d(self.survey.index,self.survey['tvdss'])
            self._northing_int = interp1d(self.survey['tvd'],self.survey.geometry.y)
            self._easting_int = interp1d(self.survey['tvd'],self.survey.geometry.x)
        else:
            self.survey=None
"""

    def to_tvd(self,md:(int,float)=None,which:list=None, ss:bool=False):
        if self._deviation is not None:
            r=[]
            _survey=self.survey
            _tvd_int = interp1d(_survey.index,_survey['tvd'])
            _tvdss_int = interp1d(_survey.index,_survey['tvdss'])

            if md is not None:
                if ss==True:
                    _tvdss = _tvdss_int(md)
                    r.append(_tvdss)
                else:
                    _tvd = _tvd_int(md)
                    r.append(_tvd)
                
            if which is not None:
                if 'perforations' in which:
                    if self._perforations is not None:
                        if ss==True:
                            self._perforations['tvdss_top']=self._perforations['md_top'].apply(_tvdss_int)
                            self._perforations['tvdss_bottom']=self._perforations['md_bottom'].apply(_tvdss_int)
                        else:
                            self._perforations['tvd_top']=self._perforations['md_top'].apply(_tvd_int)
                            self._perforations['tvd_bottom']=self._perforations['md_bottom'].apply(_tvd_int)
                            self._perforations['tvd_tick'] = self._perforations['tvd_bottom'] - self._perforations['tvd_top']
                        r.append(self._perforations)
                    else:
                        raise ValueError("No perforations have been set")


                if 'tops' in which:
                    if self._tops is not None:
                        if ss==True:
                            self._tops['tvdss_top']=self._tops['md_top'].apply(_tvdss_int)
                            self._tops['tvdss_bottom']=self._tops['md_bottom'].apply(_tvdss_int)
                        else:
                            self._tops['tvd_top']=self._tops['md_top'].apply(_tvd_int)
                            self._tops['tvd_bottom']=self._tops['md_bottom'].apply(_tvd_int)
                            self._tops['tvd_tick'] = self._tops['tvd_bottom'] - self._tops['tvd_top']
                        r.append(self._tops)
                    else:
                        raise ValueError("No tops have been set")

        else:
            raise ValueError("No survey has been set")
        return r
    

    
    def to_coord(self,md:(int,float)=None,which:list=None):
        if self._deviation is not None:  
            r=[]
            _survey=self.survey
            _northing_int = interp1d(_survey['tvd'],_survey.geometry.y)
            _easting_int = interp1d(_survey['tvd'],_survey.geometry.x)
            _tvd_int = interp1d(_survey.index,_survey['tvd'])
            if md is not None:
                _tvd = _tvd_int(md)
                _northing = _northing_int(_tvd)
                _easting = _easting_int(_tvd)
                coord = Point(_easting,_northing)
                r.append(coord)
                
            if which is not None:
                if 'perforations' in which:
                    if self._perforations is not None:
                        try:
                            self._perforations['northing'] = self._perforations['tvd_top'].apply(_northing_int)
                            self._perforations['easting'] = self._perforations['tvd_bottom'].apply(_easting_int)
                            self._perforations['geometry'] = self._perforations[['northing', 'easting']].apply(lambda x: Point(x['easting'],x['northing']),axis=1)
                            r.append(self._perforations)
                        except:
                            ValueError("No tvd has been set")
                    else:
                        raise ValueError("No perforations have been set")
                        
                if 'tops' in which:
                    if self._tops is not None:
                        try:
                            self._tops['northing'] = self._tops['tvd_top'].apply(_northing_int)
                            self._tops['easting'] = self._tops['tvd_bottom'].apply(_easting_int)
                            self._tops['geometry'] = self._tops[['northing', 'easting']].apply(lambda x: Point(x['easting'],x['northing']),axis=1)
                            r.append(self._tops)
                        except:
                            ValueError("No tvd has been set")
                    else:
                        raise ValueError("No tops have been set")
        else:
            raise ValueError("No survey has been set")
        return r

    def tops_to_logs(self,which:list=None):
        if self._tops is None:
            raise ValueError("No tops have been set")
        else:
            if which is None:
                raise ValueError("No log specification")
            else:
                if ('masterlog' in which) & (self._masterlog is not None):
                    _d = self._masterlog.df().index
                    _m = pd.DataFrame(index=_d)
                    for i in self._tops.iterrows():
                        _m.loc[(_m.index>=i[1]['md_top'])&(_m.index<=i[1]['md_bottom']),'formation'] = i[1]['formation']
                    self._masterlog.add_curve('formation',_m['formation'].values,descr='formations')
                if ('openlog' in which) & (self._openlog is not None):
                    _d = self._openlog.df().index
                    _m = pd.DataFrame(index=_d)
                    for i in self._tops.iterrows():
                        _m.loc[(_m.index>=i[1]['md_top'])&(_m.index<=i[1]['md_bottom']),'formation'] = i[1]['formation']
                    self._openlog.add_curve('formation',_m['formation'].values,descr='formations')
                if ('caselog' in which) & (self._caselog is not None):
                    _d = self._caselog.df().index
                    _m = pd.DataFrame(index=_d)
                    for i in self._tops.iterrows():
                        _m.loc[(_m.index>=i[1]['md_top'])&(_m.index<=i[1]['md_bottom']),'formation'] = i[1]['formation']
                    self._caselog.add_curve('formation',_m['formation'].values,descr='formations')

    def add_to_logs(self,df):
        if self._openlog is None:
            raise ValueError("openlog has not been added")
        else:
            #col_exist = self.openlog.df().columns
            col_add = df.columns[~df.columns.isin(np.intersect1d(df.columns, self._openlog.df().columns))]
            df_merge = self._openlog.df().merge(df[col_add], how='left', left_index=True,right_index=True)
            assert df_merge.shape[0] == self._openlog.df().shape[0]
            for i in df_merge[col_add].iteritems():
                self._openlog.add_curve(i[0],i[1])

    def interval_attributes(self,perforations:bool=False, 
                            intervals:perforations=None, 
                            curves:list = None,
                            aggfunc = ['min','max','mean']):
        if perforations == True :
            p = self._perforations
        else:
            p = intervals 
          
        curves.append('inter')
        log_appended = pd.DataFrame()
        #add column to identify the interval
        for i,c in p.iterrows():
            logdf = self._openlog.df().copy()
            logdf.loc[(logdf.index >= c['md_top'])&(logdf.index<=c['md_bottom']),'inter']=i
            
            #filter all the intervals
            logdf = logdf[~logdf['inter'].isnull()]

            #Group and aggregate functions
            log_grp = logdf[curves].groupby('inter').agg(aggfunc)
            log_appended = log_appended.append(log_grp)
        p_result = pd.concat([p,log_appended],axis=1)
        if perforations ==True:
            self._perforations = p_result 
            
        return p_result
            
        
class wells_group:
    def __init__(self,*args):
        _well_list = []

        if args is not None:
            for i in args:
                _well_list.append(i)
        
        self.wells = _well_list 

    @property
    def wells(self):
        return self._wells

    @wells.setter 
    def wells(self,value):
        assert isinstance(value,list)
        if not value:
            self._wells = {}
        else:
            assert all(isinstance(i,well) for i in value)
            w_dict={}
            for i in value:
                w_dict[i.name] = i
            self._wells = w_dict

    def add_well(self,*args):
        _add_well = []

        if args is not None:
            for i in args:
                _add_well.append(i)

        assert all(isinstance(i,well) for i in _add_well)

        _wells_dict = self.wells.copy()

        for i in _add_well:
            _wells_dict[i.name] = i
        self._wells = _wells_dict

    # Methods

    def wells_coordinates(self, wells:list=None, z_unit='ft', to_crs='EPSG:4326'):
        """
        Get a DataFrame with the wells surface coordinates
        Input:
            wells ->  (list, None) List of wells in the Group to show the matrix. 
                    If None, all wells in the group will be selected
        Return:
            wells_coord -> (pd.DataFrame) DataFrame with wells coords
        """
        assert isinstance(wells,(list,type(None)))

        # Define which wells for the distance matrix will be shown    
        if wells is None:
            _well_list = []
            for key in self.wells:
                _well_list.append(key)
        else:
            _well_list = wells

        #Create coordinates dataframe
        _coord = gpd.GeoDataFrame()

        z_coef = 0.3048 if z_unit=='ft' else 1

        for well in _well_list:
            x_coord = self.wells[well].surf_coord.x
            y_coord = self.wells[well].surf_coord.y
            z_coord = self.wells[well].surf_coord.z*z_coef if self.wells[well].surf_coord.has_z==True else 0
            shape = self.wells[well].surf_coord
            crs = self.wells[well].crs
            _w = gpd.GeoDataFrame({'x':[x_coord],'y':[y_coord],'z':[z_coord],'geometry':[shape]}, index=[well])
            _w.crs = crs
            _w = _w.to_crs(to_crs)
            _w['lon'] = _w['geometry'].x
            _w['lat'] = _w['geometry'].y
            _coord = _coord.append(_w)

        return _coord


    def wells_distance(self,wells:list=None, z=False, z_unit='ft'):
        """
        Calculate a distance matrix for the surface coordinates of the wells

        Input:
            wells ->  (list, None) List of wells in the Group to show the matrix. 
                    If None, all wells in the group will be selected
            z ->  (Bool, default False). Take into account the z component. Z must be in the same
                    units of x, y coord
            z_unit -> (str, default 'ft') Indicate the units of the z coord. 
                    If 'ft' the z is multiplied by 0.3028 otherwise by 1

        Return:
            dist_matrix -> (pd.DataFrame) Distance matrix with index and column of wells
        """
        
        
        assert isinstance(wells,(list,type(None)))

        _coord = self.wells_coordinates(wells=wells, z_unit=z_unit)

        dims = ['x','y','z'] if z==True else ['x','y']
        dist_array = distance_matrix(_coord[dims].values,_coord[dims].values)
        dist_matrix = pd.DataFrame(dist_array,index=_coord.index, columns=_coord.index)

        return dist_matrix

    def wells_map(self, wells:list=None,zoom=10, map_style = 'OpenStreetMap'):
        """
        Make a Foluim map with the selected wells

        Input:
            wells ->  (list, None) List of wells in the Group to show the matrix. 
                    If None, all wells in the group will be selected
            zoom -> (int, float) Initial zoom for folium map
        Return:
            w_map -> (folium.Map) Folium map object
        """
        assert isinstance(wells,(list,type(None)))

        _coord = self.wells_coordinates(wells=wells)

        center = _coord[['lat','lon']].mean(axis=0)

        #make the map
        map_folium = folium.Map(
            location=(center['lat'],center['lon']),
            zoom_start=zoom,
            tiles = map_style)

        for i, r in _coord.iterrows():
            folium.Marker(
                [r['lat'],r['lon']],
                tooltip=f"{i}",
                icon=folium.Icon(icon='tint', color='green')
                ).add_to(map_folium)

        return map_folium



        


        




    

