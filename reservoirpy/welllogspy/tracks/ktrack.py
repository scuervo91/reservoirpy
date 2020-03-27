import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np

def ktrack(df: pd.DataFrame, 
             k: list = None, 
             lims: list = None,
             dtick: bool =False, 
             ax=None,
             k_range:list =[0.2,20000],
             fontsize:int =8,
             correlation: pd.DataFrame = None,
             grid_numbers : list = [11,51],
             steps: list  = None,
             legend:bool = True,
             colormap: str='gist_rainbow',
             corr_kw={},
             k_kw=[]):

    #get number of curves to build the colormap
    n_curves = len(k)
    cmap = mpl.cm.get_cmap(colormap,n_curves)
    
    kax=ax or plt.gca()
    defkwa = {
    'linestyle':'-',
    'linewidth': 1
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
    if k is not None:
        for i,perm in enumerate(k):
            if len(k_kw)<i+1:
                k_kw.append(defkwa)
            k_kw[i]['color']=cmap(i)
            for (k,v) in defkwa.items():
                if k not in k_kw[i]:
                    k_kw[i][k]=v
            kax.plot(df[perm],df.index,label=perm,**k_kw[i])
    
    if lims==None: #Depth Limits
        lims=[df.index.min(),df.index.max()]

    kax.set_ylim([lims[1],lims[0]])
        
    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
     
    #Set the gridding and ticks
    kax.set_xscale("log")
    kax.set_xlim(k_range)
    ticks=np.round(np.power(10,np.linspace(np.log10(k_range[0]),np.log10(k_range[1]),int(np.log10(k_range[1]/k_range[0])+1))),decimals=1)
    kax.set_xticks(ticks)
    kax.set_xticklabels(ticks)
    kax.set_xlabel("K Perm")
    kax.xaxis.tick_top()
    kax.xaxis.set_label_position("top")
    kax.tick_params("both",labelsize=fontsize)
    kax.xformatter = 'auto'
    kax.set_yticks(minor_grid,minor=True)
    kax.set_yticks(mayor_grid)
    if dtick==True:
        kax.set_yticklabels(mayor_grid)
    else:
        kax.set_yticklabels([])
    kax.grid(True,linewidth=1.0)
    kax.grid(True,which='minor', linewidth=0.5)
    
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            kax.hlines(i[1]['depth'],k_range[0],k_range[1], **corr_kw)
            if cor_ann:
                try:
                    kax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(k_range[1]-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    kax.annotate(f"{i[1]['depth']}",xy=(k_range[1]-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
 
    if legend:
        kax.legend()