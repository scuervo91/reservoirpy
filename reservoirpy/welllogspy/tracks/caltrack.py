import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np


def caltrack(df: pd.DataFrame,
             cal: str = None, 
             bit: str = None, 
             lims: list =None,
            dtick:bool=False,
            fill:bool=False,
            fontsize: int=8,
            grid_numbers : list = [11,51],
            steps: list  = None,
            correlation: pd.DataFrame = None,
            ax=None,
            cal_kw={},
            corr_kw={},
            bit_kw={}):
    
    cal=ax or plt.gca()
    
    def_cal_kw = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 1
    }
    
    for (k,v) in def_cal_kw.items():
        if k not in cal_kw:
            cal_kw[k]=v

    def_bit_kw = {
    'color': 'darkred',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_bit_kw.items():
        if k not in bit_kw:
            bit_kw[k]=v
            
    def_corr_kw = {
    'color': 'red',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_corr_kw.items():
        if k not in corr_kw:
            corr_kw[k]=v
            
    #Set lims
    if lims==None: #Depth Limits
        lims=[df.index.min(),df.index.max()]

    cal.set_ylim([lims[1],lims[0]])
        
    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
    
    if cal is not None:
        cal.plot(df[cal],df.index,**cal_kw)
        
    if bit is not None:
        cal.plot(df[bit],df.index,**bit_kw) 
    

    cal.set_xlim([5,17])
    cal.set_xlabel("Caliper [in]")
    cal.set_xticks(np.linspace(5,17,4))
    cal.set_xticklabels(np.round(np.linspace(5,17,4),decimals=2))
    cal.xaxis.tick_top()
    cal.xaxis.set_label_position("top")
    cal.tick_params("both",labelsize=fontsize)
    cal.set_yticks(mayor_grid)
    cal.set_yticks(minor_grid,minor=True)       
    if dtick==True:
        cal.set_yticklabels(mayor_grid)
    else:
        cal.set_yticklabels([])
        
    if fill==True:
        cal.fill_betweenx(df.index,df[cal],bit,where=(cal > bit),color="orange")
        cal.fill_betweenx(df.index,df[cal],bit,where=(cal > bit),color="gray")
        
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            cal.hlines(i[1]['depth'],0,1, **corr_kw)
            if cor_ann:
                try:
                    cal.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(16-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right')
                except:
                    cal.annotate(f"{i[1]['depth']}",xy=(16-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
