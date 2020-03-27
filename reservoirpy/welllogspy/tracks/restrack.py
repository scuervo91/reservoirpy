import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np

def restrack(df: pd.DataFrame, 
             res: list = None, 
             lims: list = None,
             dtick: bool =False, 
             ax=None,
             res_range:list =[0.2,20000],
             fontsize:int =8,
             correlation: pd.DataFrame = None,
             grid_numbers : list = [11,51],
             steps: list  = None,
             legend:bool = True,
             colormap: str='summer_r',
             corr_kw={},
             res_kw=[]):
    
    #get number of curves to build the colormap
    n_curves = len(res)
    cmap = mpl.cm.get_cmap(colormap,n_curves)

    rax=ax or plt.gca()
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
    if res is not None:
        for i,r in enumerate(res):
            if len(res_kw)<i+1:
                res_kw.append(defkwa)
            res_kw[i]['color']=cmap(i)
            for (k,v) in defkwa.items():
                if k not in res_kw[i]:
                    res_kw[i][k]=v
            rax.plot(df[r],df.index,label=r,**res_kw[i])
    
    if lims==None: #Depth Limits
        lims=[df.index.min(),df.index.max()]

    rax.set_ylim([lims[1],lims[0]])
        
    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
     
    #Set the gridding and ticks
    rax.set_xscale("log")
    rax.set_xlim(res_range)
    ticks=np.round(np.power(10,np.linspace(np.log10(res_range[0]),np.log10(res_range[1]),int(np.log10(res_range[1]/res_range[0])+1))),decimals=1)
    rax.set_xticks(ticks)
    rax.set_xticklabels(ticks)
    rax.set_xlabel("Resistivity [Ohm m]")
    rax.xaxis.tick_top()
    rax.xaxis.set_label_position("top")
    rax.tick_params("both",labelsize=fontsize)
    rax.xformatter = 'auto'
    rax.set_yticks(minor_grid,minor=True)
    rax.set_yticks(mayor_grid)
    if dtick==True:
        rax.set_yticklabels(mayor_grid)
    else:
        rax.set_yticklabels([])
    rax.grid(True)
    rax.grid(True,which='minor')
    
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            rax.hlines(i[1]['depth'],res_range[0],res_range[1], **corr_kw)
            if cor_ann:
                try:
                    rax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(res_range[1]-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right')
                except:
                    rax.annotate(f"{i[1]['depth']}",xy=(res_range[1]-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
 
    if legend:
        rax.legend()