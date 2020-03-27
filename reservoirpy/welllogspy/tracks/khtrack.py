import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np


def khtrack(df:pd.DataFrame, 
            kh:str=None, 
            lims:list=None,
            dtick:bool=False,
            fill:bool=True, 
            ax=None,
            fontsize: int=8,
            grid_numbers : list = [11,51],
            steps: list  = None,
            correlation: pd.DataFrame = None,
            kh_kw={},
            corr_kw={},
            fill_kh_kw={}):
    
    hax=ax or plt.gca()
    
    def_kh_kw = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 1
    }
    
    for (k,v) in def_kh_kw.items():
        if k not in kh_kw:
            kh_kw[k]=v
            
    def_corr_kw = {
    'color': 'red',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_corr_kw.items():
        if k not in corr_kw:
            corr_kw[k]=v

    def_fill_kh_kw = {
    'color': (0.5,0.5,0.5),
    }    
    for (k,v) in def_fill_kh_kw.items():
        if k not in fill_kh_kw:
            fill_kh_kw[k]=v
    
    if kh is not None:
        hax.plot(df[kh],df.index,**kh_kw)
    
    if lims==None: #Depth Limits
        lims=[df.index.min(),df.index.max()]

    hax.set_ylim([lims[1],lims[0]])

    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
    
    hax.set_xlim([0,1])
    hax.set_xlabel("Kh Norm")
    hax.set_xticks(np.linspace(0,1,4))
    hax.set_xticklabels(np.round(np.linspace(0,1,4),decimals=2))
    hax.xaxis.tick_top()
    hax.xaxis.set_label_position("top")
    hax.tick_params("both",labelsize=fontsize)
    hax.set_yticks(mayor_grid)
    hax.set_yticks(minor_grid,minor=True)        
    if dtick==True:
        hax.set_yticklabels(mayor_grid)
    else:
        hax.set_yticklabels([])

    if fill==True:
        hax.fill_betweenx(df.index,0,kh,**fill_kh_kw)
        
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            hax.hlines(i[1]['depth'],0,1, **corr_kw)
            if cor_ann:
                try:
                    hax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right')
                except:
                    hax.annotate(f"{i[1]['depth']}",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right')
