import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np

def gastrack(df: pd.DataFrame, 
             gas: list = None, 
             lims: list = None,
             dtick: bool =False, 
             ax=None,
             gas_range:list =[0.2,20000],
             fontsize:int =8,
             correlation: pd.DataFrame = None,
             grid_numbers : list = [11,51],
             steps: list  = None,
             legend:bool = True,
             colormap: str='gist_rainbow',
             corr_kw={},
             gas_kw=[]):

    #get number of curves to build the colormap
    n_curves = len(gas)
    cmap = mpl.cm.get_cmap(colormap,n_curves)
    
    gax=ax or plt.gca()
    defkwa = {
    'linestyle':'-',
    'linewidth': 2
    }
    
    def_corr_kw = {
    'color': 'red',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_corr_kw.items():
        if k not in corr_kw:
            corr_kw[k]=v
    
    #Plot main Lines
    if gas is not None:
        for i,g in enumerate(gas):
            if len(gas_kw)<i+1:
                gas_kw.append(defkwa)
            gas_kw[i]['color']=cmap(i)
            for (k,v) in defkwa.items():
                if k not in gas_kw[i]:
                    gas_kw[i][k]=v
            gax.plot(df[g],df.index,label=g,**gas_kw[i])
    
    if lims==None: #Depth Limits
        lims=[df.index.min(),df.index.max()]

    gax.set_ylim([lims[1],lims[0]])
        
    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
     
    #Set the gridding and ticks
    gax.set_xscale("log")
    gax.set_xlim(gas_range)
    ticks=np.round(np.power(10,np.linspace(np.log10(gas_range[0]),np.log10(gas_range[1]),int(np.log10(gas_range[1]/gas_range[0])+1))),decimals=1)
    gax.set_xticks(ticks)
    gax.set_xticklabels(ticks)
    gax.set_xlabel("Gas Chromatografy")
    gax.xaxis.tick_top()
    gax.xaxis.set_label_position("top")
    gax.tick_params("both",labelsize=fontsize)
    gax.xformatter = 'auto'
    gax.set_yticks(minor_grid,minor=True)
    gax.set_yticks(mayor_grid)
    if dtick==True:
        gax.set_yticklabels(mayor_grid)
    else:
        gax.set_yticklabels([])
    gax.grid(True,linewidth=1.0)
    gax.grid(True,which='minor', linewidth=0.5)
    
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            gax.hlines(i[1]['depth'],gas_range[0],gas_range[1], **corr_kw)
            if cor_ann:
                try:
                    gax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(gas_range[1]-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right')
                except:
                    gax.annotate(f"{i[1]['depth']}",xy=(gas_range[1]-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right')
 
    if legend:
        gax.legend()