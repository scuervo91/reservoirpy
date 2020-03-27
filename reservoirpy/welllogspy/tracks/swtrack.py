import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np

def swtrack(df: pd.DataFrame,
            sw: list = None,
            lims:list =None,
            dtick: bool =False, 
            ax=None,
            fontsize:int =8,
            correlation: pd.DataFrame = None,
            grid_numbers : list = [11,51],
            steps: list  = None,
            legend:bool = True,
            fill:bool = True,
            colormap: str='winter',
            corr_kw={},
            fill_water_kw={},
            fill_oil_kw={},
            sw_kw=[]):
    
    #get number of curves to build the colormap
    n_curves = len(sw)
    cmap = mpl.cm.get_cmap(colormap,n_curves)

    sax=ax or plt.gca()
    
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

    def_fill_water_kw = {
    'color': (0.13,0.33,0.88),
    }    
    for (k,v) in def_fill_water_kw.items():
        if k not in fill_water_kw:
            fill_water_kw[k]=v

    def_fill_oil_kw = {
    'color': (0.07,0.72,0.08),
    }    
    for (k,v) in def_fill_oil_kw.items():
        if k not in fill_oil_kw:
            fill_oil_kw[k]=v

    #Plot main Lines
    if sw is not None:
        for i,r in enumerate(sw):
            if len(sw_kw)<i+1:
                sw_kw.append(defkwa)
            sw_kw[i]['color']=cmap(i)
            for (k,v) in defkwa.items():
                if k not in sw_kw[i]:
                    sw_kw[i][k]=v
            sax.plot(df[r],df.index,label=r,**sw_kw[i])
    
    if lims==None: #Depth Limits
        lims=[df.index.min(),df.index.max()]

    sax.set_ylim([lims[1],lims[0]])
        
    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
    
    sax.set_xlim([1,0])
    sax.set_xlabel("Sw")
    sax.set_xticks(np.linspace(0,1,4))
    sax.set_xticklabels(np.round(np.linspace(0,1,4),decimals=2))
    sax.xaxis.tick_top()
    sax.xaxis.set_label_position("top")
    sax.tick_params("both",labelsize=fontsize)
    sax.set_yticks(minor_grid,minor=True)
    sax.set_yticks(mayor_grid)
    if dtick==True:
        sax.set_yticklabels(mayor_grid)
    else:
        sax.set_yticklabels([])
    if fill==True:
        sax.fill_betweenx(df.index,1,df[sw[0]],**fill_oil_kw)
        sax.fill_betweenx(df.index,df[sw[0]],0,**fill_water_kw)
        
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            sax.hlines(i[1]['depth'],0,1, **corr_kw)
            if cor_ann:
                try:
                    sax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right')
                except:
                    sax.annotate(f"{i[1]['depth']}",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right')



