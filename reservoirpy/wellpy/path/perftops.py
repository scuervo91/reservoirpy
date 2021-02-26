import geopandas as gpd 


class Perforations(gpd.GeoDataFrame):
    """Perforations [GeoDataFrame subclass that represent the perforations  
                    intervals if a well. In order an instance of this class
                    be functional, a certain columns names must be intdicating when existing.
                    
                    To refer a perforations top and bottom in Measured Depth (MD) the columns are:
                        - 'md_top' 'md_bottom'

                    To refer a perforations top and bottom in True Vertical Depth (TVD) the columns are:
                        - 'tvd_top' 'tvd_bottom'
                        
                    To refer a perforations top and bottom in True Vertical Depth Subsea (TVDSS) the columns are:
                        - 'tvdss_top' 'tvdss_bottom'
                        
                    As It is a GeoDataFrame subclass, it can be added a 'geometry' column with the 
                    shapely data. https://geopandas.org/
                    
                    ]

    Methods
    ----------
    get_tick : 
        [If present the right columns names, it tries to estimate the thickness in md and tvd
        
        If successful the return columns are 'md_tick' and 'tvd_tick'        
        ]
        
    get_mid_point : 
        [If present the right columns names, it tries to estimate the mid point in md and tvd
        
        If successful the return columns are 'md_mid_point' and 'tvd_mid_point'        
        ]
    """
    def __init__(self, *args, **kwargs):
        kh = kwargs.pop("kh", None)
        productivity_index = kwargs.pop("productivity_index", None)
        is_open = kwargs.pop('is_open',None)  
        fluid = kwargs.pop('fluid',None)                                                                                                                               
        super(Perforations, self).__init__(*args, **kwargs)

        if kh is not None:
            assert isinstance(kh,list) 
            assert all(isinstance(i,(int,float)) for i in kh)
            self['kh'] = kh
        elif 'kh' in self.columns:
            assert all(isinstance(i,(int,float)) for i in self['kh'].tolist())

        if productivity_index is not None:
            assert isinstance(productivity_index,list) 
            assert all(isinstance(i,(int,float)) for i in productivity_index)
            self['productivity_index'] = productivity_index
        elif 'productivity_index' in self.columns:
            assert all(isinstance(i,(int,float)) for i in self['productivity_index'].tolist())

        if is_open is not None:
            assert isinstance(is_open,list)
            assert all(isinstance(i,bool) for i in is_open)
            self['is_open'] = is_open
        elif 'is_open' in self.columns:
            assert all(isinstance(i,bool) for i in self['is_open'].tolist())

        if fluid is not None:
            assert isinstance(fluid,list)
            assert all(i in ['oil','gas','water'] for i in fluid)
            self['fluid'] = fluid
        elif 'fluid' in self.columns:
            assert all(isinstance(i,bool) for i in self['is_open'].tolist())
            assert all(i in ['oil','gas','water'] for i in self['fluid'].tolist())


    def open_perf(self):
        return self[self['is_open']==True]

    def get_tick(self):
        """get_tick [Estimate the tick in md and tvd.]

        Returns
        -------
        [Perforations]
            [Perforations GeoDataFrame with the corresponding tick columns]
        """
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
        """get_mid_point [Estimate the tick in md and tvd.]

        Returns
        -------
        [Perforations]
            [Perforations GeoDataFrame with the corresponding mid_point columns]
        """
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
        return Perforations
    
class Tops(gpd.GeoDataFrame):
    """Tops [GeoDataFrame subclass that represent the well formation tops 
                    intervals if a well. In order an instance of this class
                    be functional, a certain columns names must be intdicating when existing.
                    
                    To refer a well formation top and bottom in Measured Depth (MD) the columns are:
                        - 'md_top' 'md_bottom'

                    To refer a well formation top and bottom in True Vertical Depth (TVD) the columns are:
                        - 'tvd_top' 'tvd_bottom'
                        
                    To refer a well formation top and bottom in True Vertical Depth Subsea (TVDSS) the columns are:
                        - 'tvdss_top' 'tvdss_bottom'
                        
                    As It is a GeoDataFrame subclass, it can be added a 'geometry' column with the 
                    shapely data. https://geopandas.org/
                    
                    ]

    Methods
    ----------
    get_tick : 
        [If present the right columns names, it tries to estimate the thickness in md and tvd
        
        If successful the return columns are 'md_tick' and 'tvd_tick'        
        ]
        
    get_mid_point : 
        [If present the right columns names, it tries to estimate the mid point in md and tvd
        
        If successful the return columns are 'md_mid_point' and 'tvd_mid_point'        
        ]
    """
    def __init__(self, *args, **kwargs):    
        formation = kwargs.pop("formation", None)            
        unit = kwargs.pop("unit", None)                                                                                                                     
        super(Tops, self).__init__(*args, **kwargs)  

        if formation is not None:
            formation = np.atleast_1d(formation)
            self['formation'] = formation
            self.set_index('formation',inplace=True)
        elif 'formation' in self.columns:
            self.set_index('formation',inplace=True)

        if unit is not None:
            unit = np.atleast_1d(unit)
            self['unit'] = unit
            self.set_index('unit',inplace=True)
        elif 'unit' in self.columns:
            self.set_index('unit',inplace=True)

    def get_tick(self):
        """get_tick [Estimate the tick in md and tvd.]

        Returns
        -------
        [Perforations]
            [Perforations GeoDataFrame with the corresponding tick columns]
        """
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
        """get_mid_point [Estimate the tick in md and tvd.]

        Returns
        -------
        [Perforations]
            [Perforations GeoDataFrame with the corresponding mid_point columns]
        """
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
        return Tops