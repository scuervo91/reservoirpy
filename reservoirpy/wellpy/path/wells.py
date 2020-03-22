import pandas as pd 
import numpy as np 
from .checkarrays import checkarrays, checkarrays_tvd, checkarrays_monotonic_tvd
from .mincurve import minimum_curvature, min_curve_method
from .interpolate import interpolate_deviation, interpolate_position
from scipy.interpolate import interp1d
from shapely.geometry import Point

class perforation:
    def __init__(self,top,base,depth_type,skin,rw):
        self.top = top
        self.base = base 
        self.depth_type = depth_type
        self.skin = skin 
        self.rw = rw 

class well:
    def __init__(self,  *args, **kwargs):
        self.name = kwargs.pop("name", None)
        self.rte = kwargs.pop("rte", None)
        self.surf_coord = kwargs.pop("surf_coord", None)
        self.perforations = kwargs.pop("perforations",None)
        _deviation = pd.DataFrame(kwargs.pop("deviation",None))
        self.survey = min_curve_method(_deviation['md'],_deviation['inc'],_deviation['azi'],surface_easting=self.surf_coord.x, surface_northing=self.surf_coord.y, kbe=self.rte)
        
    def sample_deviation(self,step=100):
        new_dev = interpolate_deviation(self.survey.index, self.survey['inc'], self.survey['azi'], md_step=step)
        return new_dev

    def sample_position(self,step=100):
        new_pos = interpolate_position(self.survey['tvd'], self.survey['easting'], self.survey['northing'], tvd_step=step)
        return new_pos
    
    def to_tvd(self,md):
        tvd_int = interp1d(self.survey.index,self.survey['tvd'])
        tvd = tvd_int(md)
        return tvd 
    
    def to_tvdss(self,md):
        tvdss_int = interp1d(self.survey.index,self.survey['tvdss'])
        tvdss = tvdss_int(md)
        return tvdss
