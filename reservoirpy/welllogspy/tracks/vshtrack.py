import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np

def vshtrack(df: pd.DataFrame, 
            vsh: str = None, 
            lims: list = None,
            inverse: bool = False,
            dtick: bool = False,
            fill: bool = True,
            fontsize: int = 8,
            correlation: pd.DataFrame = None,
            grid_numbers : list = [11,51],
            steps: list  = None,
            vsh_kw={},
            corr_kw={},
            ax=None
            ):
    """
    

    Parameters
    ----------
    df : pd.DataFrame
        DESCRIPTION.
    vsh : str, optional
        DESCRIPTION. The default is None.
    lims : list, optional
        DESCRIPTION. The default is None.
    inverse : bool, optional
        DESCRIPTION. The default is False.
    dtick : bool, optional
        DESCRIPTION. The default is False.
    fill : bool, optional
        DESCRIPTION. The default is False.
    cutoff : float, optional
        DESCRIPTION. The default is 0.5.
    fontsize : int, optional
        DESCRIPTION. The default is 8.
    grid_numbers : list, optional
        DESCRIPTION. The default is [11,51].
    steps : list, optional
        DESCRIPTION. The default is None.
    vsh_kw : TYPE, optional
        DESCRIPTION. The default is {}.
    ax : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """

    
    vax=ax or plt.gca()
    def_vsh_kw = {
    'color': 'black',
    'linestyle':'--',
    'linewidth': 1
    }
    
    for (k,v) in def_vsh_kw.items():
        if k not in vsh_kw:
            vsh_kw[k]=v

    def_corr_kw = {
    'color': 'red',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_corr_kw.items():
        if k not in corr_kw:
            corr_kw[k]=v
    
    if vsh is not None:
        vax.plot(df[vsh],df.index,**vsh_kw)
        
    # Set Axis parameters 
    if lims==None: #Depth Limits
        lims=[df.index.min(),df.index.max()]
    
    vax.set_ylim([lims[1],lims[0]])
        
    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
    
    #Set the gridding and ticks
    vax.set_xlim([1,0]) if inverse==True else vax.set_xlim([0,1])
    vax.set_xlabel("Vshale")
    vax.set_xticks(np.linspace(0,1,4))
    vax.set_xticklabels(np.round(np.linspace(0,1,4),decimals=2))
    vax.xaxis.tick_top()
    vax.xaxis.set_label_position("top")
    vax.tick_params("both",labelsize=fontsize)
    vax.set_yticks(minor_grid,minor=True)
    vax.set_yticks(mayor_grid)
    if dtick==True:
        vax.set_yticklabels(mayor_grid)
    else:
        vax.set_yticklabels([])
        
    if fill==True:
       # cmap=mpl.cm.get_cmap('YlOrBr',128)
        cmap=mpl.colors.LinearSegmentedColormap.from_list('vsh',
          [(252/256,210/256,23/256),(66/256,59/256,20/256),(41/256,29/256,13/256)],N=256)
        df_filtered = df.loc[(df.index>=lims[0])&(df.index<=lims[1]),vsh]
        for i in range(0,df_filtered.shape[0]-1):
            vax.fill_betweenx([df_filtered.index[i],df_filtered.index[i+1]],
                              [df_filtered.iloc[i],df_filtered.iloc[i+1]],
                              1,color=cmap(df_filtered.iloc[i+1]))
    
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            vax.hlines(i[1]['depth'],0,1, **corr_kw)
            if cor_ann:
                try:
                    vax.annotate(f"{i[1]['depth']} - {i[1]['comment']}",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right')
                except:
                    vax.annotate(f"{i[1]['depth']}",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
     
    
    if 'label' in vsh_kw:
        vax.legend()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        