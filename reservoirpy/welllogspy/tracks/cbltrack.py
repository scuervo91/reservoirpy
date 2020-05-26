import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np


def cbltrack(df:pd.DataFrame,
             cbl: str=None, 
             lims: list=None,
             dtick: bool=False, 
             ax=None,
             fontsize=8,
             correlation: pd.DataFrame = None,
             grid_numbers : list = [11,51],
             steps: list  = None,
             legend:bool = True,
             cbl_colormap:str='cool',
             cbl_kw:dict={},
             cblx5_kw:dict={},
             corr_kw:dict={},
             depth_ref:str='md'):
    
    cblax=ax or plt.gca()
    cbl2ax=cblax.twiny()
    
    def_cbl_kw = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 1
    }
    
    for (k,v) in def_cbl_kw.items():
        if k not in cbl_kw:
            cbl_kw[k]=v

    def_cblx5_kw = {
    'color': 'blue',
    'linestyle':'-',
    'linewidth': 2
    }
    
    for (k,v) in def_cblx5_kw.items():
        if k not in cblx5_kw:
            cblx5_kw[k]=v

    depth = df.index if depth_ref=='md' else df[depth_ref]
    if cbl is not None:
        cblax.plot(df[cbl],df.index,**cbl_kw)
        cbl2ax.plot(df[cbl],df.index,**cblx5_kw)
        if isinstance(cbl,str):
            cblax.plot(df[cbl],depth,**cbl_kw)   #Plotting
            cbl2ax.plot(df[cbl],depth,**cblx5_kw)
        elif isinstance(cbl,list):
            cmap = mpl.cm.get_cmap(cbl_colormap,len(cbl))
            for i,g in enumerate(cbl):
                cbl_kw['color']=cmap(i)
                cblax.plot(df[g],depth,**cbl_kw)
                cbl2ax.plot(df[g],depth,**cblx5_kw)
    
    # Set The lims of depth               
    if lims==None: #Depth Limits
        lims=[depth.min(),depth.max()]

    cblax.set_ylim([lims[1],lims[0]])
        
    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
    
    cblax.set_xlim([0,100])
    cbl2ax.set_xlim([0,20])
    cblax.set_xlabel("CBL")
    cbl2ax.set_xlabel("CBLx5")
    cblax.set_xticks(np.linspace(0,100,4))
    cbl2ax.set_xticks(np.linspace(0,20,4))
    cblax.set_xticklabels(np.round(np.linspace(0,100,4),decimals=2))
    cbl2ax.set_xticklabels(np.round(np.linspace(0,20,4),decimals=2))
    cblax.xaxis.tick_top()
    cblax.xaxis.set_label_position("top")
    cblax.tick_params("both",labelsize=fontsize)
    cbl2ax.tick_params("both",labelsize=fontsize)
    cblax.set_yticks(mayor_grid)
    cblax.set_yticks(minor_grid,minor=True)        
    cblax.grid(True,linewidth=1.0)
    cblax.grid(True,which='minor', linewidth=0.5)
    if dtick==True:
        cblax.set_yticklabels(mayor_grid)
    else:
        cblax.set_yticklabels([])
        
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            cblax.hlines(i[1]['depth'],0,100, **corr_kw)
            if cor_ann:
                try:
                    cblax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(100-30,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    cblax.annotate(f"{i[1]['depth']}",xy=(100-30,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})