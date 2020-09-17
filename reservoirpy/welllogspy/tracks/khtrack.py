import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np


def khtrack(df:pd.DataFrame, 
            kh:(list,str)=None, 
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
            fill_kh_kw={},
            kh_colormap:str='gray',
            depth_ref:str='md'):
    """khtrack [summary]

    Parameters
    ----------
    df : pd.DataFrame
        [description]
    kh : [type], optional
        [description], by default None
    lims : list, optional
        [description], by default None
    dtick : bool, optional
        [description], by default False
    fill : bool, optional
        [description], by default True
    ax : [type], optional
        [description], by default None
    fontsize : int, optional
        [description], by default 8
    grid_numbers : list, optional
        [description], by default [11,51]
    steps : list, optional
        [description], by default None
    correlation : pd.DataFrame, optional
        [description], by default None
    kh_kw : dict, optional
        [description], by default {}
    corr_kw : dict, optional
        [description], by default {}
    fill_kh_kw : dict, optional
        [description], by default {}
    kh_colormap : str, optional
        [description], by default 'gray'
    depth_ref : str, optional
        [description], by default 'md'
    """
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

    depth = df.index if depth_ref=='md' else df[depth_ref]
    if kh is not None:
        if isinstance(kh,str):
            hax.plot(df[kh],depth,**kh_kw)   #Plotting
        elif isinstance(kh,list):
            cmap = mpl.cm.get_cmap(kh_colormap,len(kh))
            for i,g in enumerate(kh):
                kh_kw['color']=cmap(i)
                hax.plot(df[g],depth,**kh_kw)
    
    if lims==None: #Depth Limits
        lims=[depth.min(),depth.max()]

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
        hax.fill_betweenx(depth,0,df[kh],**fill_kh_kw)
        
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            hax.hlines(i[1]['depth'],0,1, **corr_kw)
            if cor_ann:
                try:
                    hax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    hax.annotate(f"{i[1]['depth']}",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
