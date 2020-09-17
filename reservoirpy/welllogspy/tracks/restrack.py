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
             corr_kw:dict={},
             res_kw:dict=[],
             depth_ref:str='md'):
    """restrack [summary]

    Parameters
    ----------
    df : pd.DataFrame
        [description]
    res : list, optional
        [description], by default None
    lims : list, optional
        [description], by default None
    dtick : bool, optional
        [description], by default False
    ax : [type], optional
        [description], by default None
    res_range : list, optional
        [description], by default [0.2,20000]
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
        [description], by default 'summer_r'
    corr_kw : dict, optional
        [description], by default {}
    res_kw : dict, optional
        [description], by default []
    depth_ref : str, optional
        [description], by default 'md'
    """
    assert isinstance(df,pd.DataFrame)
    assert depth_ref in ['md','tvd','tvdss'], "depth_ref can only be one of ['md','tvd','tvdss']"
    
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

    depth = df.index if depth_ref=='md' else df[depth_ref]
    #Plot main Lines
    if res is not None:
        for i,r in enumerate(res):
            if len(res_kw)<i+1:
                res_kw.append(defkwa)
            res_kw[i]['color']=cmap(i)
            for (k,v) in defkwa.items():
                if k not in res_kw[i]:
                    res_kw[i][k]=v
            rax.plot(df[r],depth,label=r,**res_kw[i])
    
    if lims==None: #Depth Limits
        lims=[depth.min(),depth.max()]

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
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    rax.annotate(f"{i[1]['depth']}",xy=(res_range[1]-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
 
    if legend:
        rax.legend()