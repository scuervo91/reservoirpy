import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np

def oilshowtrack(df: pd.DataFrame,
                 oilshow: str = None, 
                 lims: list = None,
                 dtick: bool = False,
                 fill: bool = True, 
                 ax=None,
                 fontsize=8,
                 correlation: pd.DataFrame = None,
                 grid_numbers : list = [11,51],
                 steps: list  = None,
                 corr_kw={},
                 show_kw={},
                 fill_kw={}):

    oax=ax or plt.gca()
    
    defkwa = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 1
    }
    
    for (k,v) in defkwa.items():
        if k not in show_kw:
            show_kw[k]=v
    
    def_corr_kw = {
    'color': 'red',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_corr_kw.items():
        if k not in corr_kw:
            corr_kw[k]=v

    def_fill_kw = {
    'color': 'darkgreen',
    }    
    for (k,v) in def_fill_kw.items():
        if k not in fill_kw:
            fill_kw[k]=v
            
    if oilshow is not None:
        oax.plot(df[oilshow],df.index,**show_kw)
    
    if lims==None: #Depth Limits
        lims=[df.index.min(),df.index.max()]

    oax.set_ylim([lims[1],lims[0]])

    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
    
    oax.set_xlim([0,1])
    oax.set_xlabel("OilShow")
    oax.set_xticks(np.linspace(0,1,4))
    oax.set_xticklabels(np.round(np.linspace(0,1,4),decimals=2))
    oax.xaxis.tick_top()
    oax.xaxis.set_label_position("top")
    oax.tick_params("both",labelsize=fontsize)
    oax.set_yticks(mayor_grid)
    oax.set_yticks(minor_grid,minor=True)        
    if dtick==True:
        oax.set_yticklabels(mayor_grid,11)
    else:
        oax.set_yticklabels([])
    if fill==True:
        oax.fill_betweenx(df.index,0,df[oilshow],**fill_kw)
        
        
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            oax.hlines(i[1]['depth'],0,1, **corr_kw)
            if cor_ann:
                try:
                    oax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    oax.annotate(f"{i[1]['depth']}",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
    