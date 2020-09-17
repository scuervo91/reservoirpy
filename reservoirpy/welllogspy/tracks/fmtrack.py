import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np

def fmtrack(df,
            lims:list = None,
            depth_ref:str = 'md',
            area_kw:dict = {},
            line_kw:dict = {},
            fm_colormap: str='winter',
            ax=None
            ):
    """fmtrack [summary]

    Parameters
    ----------
    df : [type]
        [description]
    lims : list, optional
        [description], by default None
    depth_ref : str, optional
        [description], by default 'md'
    area_kw : dict, optional
        [description], by default {}
    line_kw : dict, optional
        [description], by default {}
    fm_colormap : str, optional
        [description], by default 'winter'
    ax : [type], optional
        [description], by default None
    """
    
    #fm=None,lims=None, top='top',bottom='bottom',name='name',ax=None):

    # Default kwargs for area
    def_area_kw = {
    'xmin': 0,
    'xmax': 1,
    'alpha': 0.2
    }
    for (k,v) in def_area_kw.items():
        if k not in area_kw:
            area_kw[k]=v

    # Default kwargs for lines
    def_line_kw = {
        'linewidth': 1
    }
    for (k,v) in def_line_kw.items():
        if k not in line_kw:
            line_kw[k]=v

    fax=ax or plt.gca()
        
    if lims==None: #Depth Limits
        lims=[df.loc[:,f'{depth_ref}_top'].min(),df.loc[:,f'{depth_ref}_bottom'].max()]
    fax.set_ylim([lims[1],lims[0]])

    dff = df[(df[f'{depth_ref}_top']<=lims[1]) & (df[f'{depth_ref}_bottom']>=lims[0])]
    
    for i in dff.iterrows():
        fax.axhspan(i[1][f'{depth_ref}_top'],i[1][f'{depth_ref}_bottom'],**area_kw)
        fax.hlines([i[1][f'{depth_ref}_top'],i[1][f'{depth_ref}_bottom']],0,1, line_kw)
        
    fax.set_yticks(dff[f'{depth_ref}_top'])
    fax.set_yticklabels(dff['formation'])