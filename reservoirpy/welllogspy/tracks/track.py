import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np

def track(df: pd.DataFrame,
          curves: list = None,
          lims: list = None, 
          clims: list = None,
          dtick: bool = False,
          scale: str ='linear',
          curvetitle: str ='Track', 
          colormap: str='plasma',
          ax=None,
          fontsize=8,
          correlation: pd.DataFrame = None,
          grid_numbers : list = [11,51],
          steps: list  = None,
          legend:bool = True,
          grid:bool = True,
          track_kw=[],
          corr_kw={},
          depth_ref:str='md',
          ):
    """track [summary]

    Parameters
    ----------
    df : pd.DataFrame
        [description]
    curves : list, optional
        [description], by default None
    lims : list, optional
        [description], by default None
    clims : list, optional
        [description], by default None
    dtick : bool, optional
        [description], by default False
    scale : str, optional
        [description], by default 'linear'
    curvetitle : str, optional
        [description], by default 'Track'
    colormap : str, optional
        [description], by default 'plasma'
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
    grid : bool, optional
        [description], by default True
    track_kw : list, optional
        [description], by default []
    corr_kw : dict, optional
        [description], by default {}
    depth_ref : str, optional
        [description], by default 'md'
    """
    #get number of curves to build the colormap
    n_curves = len(curves)
    cmap = mpl.cm.get_cmap(colormap,n_curves)

    tax=ax or plt.gca()
    
    defkwa = {
    'color': 'black',
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
    if curves is not None:
        for i,g in enumerate(curves):
            if len(track_kw)<i+1:
                track_kw.append(defkwa)
                track_kw[i]['color']=cmap(i)
            for (k,v) in defkwa.items():
                if k not in track_kw[i]:
                    track_kw[i][k]=v

            scatter = track_kw[i].pop('scatter',False)
            if scatter:
                tax.scatter(df[g],depth,label=g,**track_kw[i])
            else:
                tax.plot(df[g],depth,label=g,**track_kw[i])
    
    if lims==None: #Depth Limits
        lims=[depth.min(),depth.max()]

    tax.set_ylim([lims[1],lims[0]])

    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])


    tax.set_xlabel(curvetitle)
    tax.set_xscale(scale)
    if clims==None: 
        clims=[0,1]

    tax.set_xlim(clims)

    if scale=='log':
        ticks=np.round(np.power(10,np.linspace(np.log10(clims[0]),np.log10(clims[1]),int(np.log10(clims[1]/clims[0])+1))),decimals=1)
    else:
        ticks = np.round(np.linspace(clims[0],clims[1],4),decimals=1)

    tax.set_xticks(ticks)
    tax.set_xticklabels(ticks)
    tax.xaxis.tick_top()
    tax.xaxis.set_label_position("top")
    tax.tick_params("both",labelsize=fontsize)
    tax.set_yticks(mayor_grid)
    tax.set_yticks(minor_grid,minor=True)    

    if grid == True:
        tax.grid(True,linewidth=1.0)
        tax.grid(True,which='minor', linewidth=0.5)
        
    if dtick==True:
        tax.set_yticklabels(np.linspace(lims[0],lims[1],11))
    else:
        tax.set_yticklabels([])
    if legend:
        tax.legend() 
        
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            tax.hlines(i[1]['depth'],clims[0],clims[1], **corr_kw)
            if cor_ann:
                try:
                    tax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(clims[1]-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    tax.annotate(f"{i[1]['depth']}",xy=(clims[1]-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
 
    