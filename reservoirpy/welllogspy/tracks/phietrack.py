import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np

def phietrack(df: pd.DataFrame,
              phi: list =None, 
             lims: list = None,
             phi_range :list = [0,0.35],
             dtick: bool =False, 
             ax=None,
             fontsize=8,
             correlation: pd.DataFrame = None,
             grid_numbers : list = [11,51],
             steps: list  = None,
             legend:bool = True,
             colormap: str='Dark2',
             corr_kw={},
             phi_kw:list = [],
             depth_ref:str='md'):
    """phietrack [summary]

    Parameters
    ----------
    df : pd.DataFrame
        [description]
    phi : list, optional
        [description], by default None
    lims : list, optional
        [description], by default None
    phi_range : list, optional
        [description], by default [0,0.35]
    dtick : bool, optional
        [description], by default False
    ax : [type], optional
        [description], by default None
    fontsize : int, optional
        [description], by default 8
    correlation : pd.DataFrame, optional
        [description], by default None
    grid_numbers : list, optional
        [description], by default [11,51]
    steps : list, optional
        [description], by default None
    legend : bool, optional
        [description], by default True
    colormap : str, optional
        [description], by default 'Dark2'
    corr_kw : dict, optional
        [description], by default {}
    phi_kw : list, optional
        [description], by default []
    depth_ref : str, optional
        [description], by default 'md'
    """
    #get number of curves to build the colormap
    n_curves = len(phi)
    cmap = mpl.cm.get_cmap(colormap,n_curves)
    
    pax=ax or plt.gca()
    
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

    depth = df.index if depth_ref=='md' else df[depth_ref]
    #Plot main Lines
    if phi is not None:
        for i,r in enumerate(phi):
            if len(phi_kw)<i+1:
                phi_kw.append(defkwa)
            phi_kw[i]['color']=cmap(i)
            for (k,v) in defkwa.items():
                if k not in phi_kw[i]:
                    phi_kw[i][k]=v
            pax.plot(df[r],depth,label=r,**phi_kw[i])
    
    if lims==None: #Depth Limits
        lims=[depth.min(),depth.max()]

    pax.set_ylim([lims[1],lims[0]])

    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
    
    pax.set_xlim(phi_range)
    pax.set_xlabel("phie")
    pax.set_xticks(np.linspace(phi_range[0],phi_range[1],4))
    pax.set_xticklabels(np.round(np.linspace(phi_range[0],phi_range[1],4),decimals=2))
    pax.xaxis.tick_top()
    pax.xaxis.set_label_position("top")
    pax.tick_params("both",labelsize=fontsize)
    pax.set_yticks(minor_grid,minor=True)
    pax.set_yticks(mayor_grid)       
    pax.grid(True,linewidth=1.0)
    pax.grid(True,which='minor', linewidth=0.5)
    if dtick==True:
        pax.set_yticklabels(mayor_grid)
    else:
        pax.set_yticklabels([])

    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            pax.hlines(i[1]['depth'],0,1, **corr_kw)
            if cor_ann:
                try:
                    pax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(0.35-0.05,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    pax.annotate(f"{i[1]['depth']}",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})

        
    if legend:
        pax.legend()