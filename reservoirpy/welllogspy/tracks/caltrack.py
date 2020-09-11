import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np


def caltrack(df: pd.DataFrame,
             cali: (list,str) = None, 
             bit: (list,str) = None, 
             lims: (list) =None,
             cal_lim:list=[5,20],
            dtick:bool=False,
            fill:bool=False,
            fontsize: int=8,
            grid_numbers : list = [11,51],
            steps: list  = None,
            correlation: pd.DataFrame = None,
            ax=None,
            cal_kw:dict={},
            corr_kw:dict={},
            bit_kw:dict={},
            depth_ref:str='md',
            cal_colormap: str='winter',
            bit_colormap: str='bone'):

    assert isinstance(df,pd.DataFrame)
    assert depth_ref in ['md','tvd','tvdss'], "depth_ref can only be one of ['md','tvd','tvdss']"

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
    
    depth = df.index if depth_ref=='md' else df[depth_ref] 

    if cali is not None:
        if isinstance(cali,str):
            cal.plot(df[cali],depth,**cal_kw)   #Plotting
        elif isinstance(cali,list):
            cmap = mpl.cm.get_cmap(cal_colormap,len(cal))
            for i,c in enumerate(cal):
                cal_kw['color']=cmap(i)
                cal.plot(df[c],depth,**cal_kw)
        
    if bit is not None:
        cal.plot(df[bit],depth,**bit_kw) 
    

    cal.set_xlim(cal_lim)
    cal.set_xlabel("Caliper [in]")
    cal.set_xticks(np.linspace(cal_lim[0],cal_lim[1],4))
    cal.set_xticklabels(np.round(np.linspace(cal_lim[0],cal_lim[1],4),decimals=1))
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
        cal.fill_betweenx(depth,df[cali],df[bit],where=(cali > bit),color="orange")
        cal.fill_betweenx(depth,df[cali],df[bit],where=(cali < bit),color="gray")
        
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            cal.hlines(i[1]['depth'],0,1, **corr_kw)
            if cor_ann:
                try:
                    cal.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(16-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    cal.annotate(f"{i[1]['depth']}",xy=(16-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
