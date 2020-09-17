import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np

def vshtrack(df: pd.DataFrame, 
            vsh: (list,str) = None, 
            lims: list = None,
            inverse: bool = False,
            dtick: bool = False,
            fill: bool = True,
            fontsize: int = 8,
            correlation: pd.DataFrame = None,
            fm: pd.DataFrame = None,
            perf: pd.DataFrame = None,
            grid_numbers : list = [11,51],
            steps: list  = None,
            vsh_kw:dict={},
            corr_kw:dict={},
            fm_kw:dict={},
            perf_kw:dict={},
            vsh_colormap:str='copper',
            depth_ref:str='md',
            ax=None
            ):
    """vshtrack [summary]

    Parameters
    ----------
    df : pd.DataFrame
        [description]
    vsh : [type], optional
        [description], by default None
    lims : list, optional
        [description], by default None
    inverse : bool, optional
        [description], by default False
    dtick : bool, optional
        [description], by default False
    fill : bool, optional
        [description], by default True
    fontsize : int, optional
        [description], by default 8
    correlation : pd.DataFrame, optional
        [description], by default None
    fm : pd.DataFrame, optional
        [description], by default None
    perf : pd.DataFrame, optional
        [description], by default None
    grid_numbers : list, optional
        [description], by default [11,51]
    steps : list, optional
        [description], by default None
    vsh_kw : dict, optional
        [description], by default {}
    corr_kw : dict, optional
        [description], by default {}
    fm_kw : dict, optional
        [description], by default {}
    perf_kw : dict, optional
        [description], by default {}
    vsh_colormap : str, optional
        [description], by default 'copper'
    depth_ref : str, optional
        [description], by default 'md'
    ax : [type], optional
        [description], by default None
    """
    assert isinstance(df,pd.DataFrame)
    assert depth_ref in ['md','tvd','tvdss'], "depth_ref can only be one of ['md','tvd','tvdss']"
    
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

    def_fm_kw = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 2
    }    
    for (k,v) in def_fm_kw.items():
        if k not in fm_kw:
            fm_kw[k]=v

    def_perf_kw = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 1
    }    
    for (k,v) in def_perf_kw.items():
        if k not in perf_kw:
            perf_kw[k]=v

    depth = df.index if depth_ref=='md' else df[depth_ref]
    if vsh is not None:
        if isinstance(vsh,str):
            vax.plot(df[vsh],depth,**vsh_kw)   #Plotting
        elif isinstance(vsh,list):
            cmap = mpl.cm.get_cmap(vsh_colormap,len(vsh))
            for i,g in enumerate(vsh):
                vsh_kw['color']=cmap(i)
                vax.plot(df[g],depth,**vsh_kw)
        
    # Set Axis parameters 
    if lims==None: #Depth Limits
        lims=[df.depth.min(),df.depth.max()]
    
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
        
    if fill==True and isinstance(vsh,str):
       # cmap=mpl.cm.get_cmap('YlOrBr',128)
        cmap=mpl.colors.LinearSegmentedColormap.from_list('vsh',
          [(252/256,210/256,23/256),(66/256,59/256,20/256),(41/256,29/256,13/256)],N=256)
        df_filtered = df.loc[(depth>=lims[0])&(depth<=lims[1]),vsh]
        depth_f = df_filtered.index if depth_ref=='md' else df_filtered[depth_ref]
        for i in range(0,df_filtered.shape[0]-1):
            vax.fill_betweenx([depth_f[i],depth_f[i+1]],
                              [df_filtered.iloc[i],df_filtered.iloc[i+1]],
                              1,color=cmap(df_filtered.iloc[i+1]))

    #Add Formations tops
    if fm is not None:
        fm_ann = fm_kw.pop('ann',False)
        for i in fm.iterrows():
            if i[1][f'{depth_ref}_top'] < lims[0] or i[1][f'{depth_ref}_top'] > lims[1]:
                continue
            vax.hlines([i[1][f'{depth_ref}_top']],0,1, **fm_kw)
            if fm_ann:
               vax.annotate(f"Top of {i[0]}",xy=(1-0.7,i[1][f'{depth_ref}_top']-2),
                             xycoords='data',horizontalalignment='right',
                             bbox={'boxstyle':'round', 'fc':'0.8'})
                          

     #Add Interval Perforated
    if perf is not None:
        perf_ann = perf_kw.pop('ann',False)
        for i in perf.iterrows():
            if perf_ann:
                try:
                    vax.annotate(f"Top:{i[1][f'{depth_ref}_top']} \nBottom:{i[1][f'{depth_ref}_bottom']} \nNote:{i[1]['comment']}",
                                  xy=(0,(i[1][f'{depth_ref}_top']+i[1][f'{depth_ref}_bottom'])/2),xycoords='data',
                                  xytext=(-180, 0), textcoords='offset points',
                                  arrowprops=dict(arrowstyle="->"))
                except:
                    vax.annotate(f"Top:{i[1][f'{depth_ref}_top']} \nBottom:{i[1][f'{depth_ref}_bottom']}",
                    xy=(0,(i[1][f'{depth_ref}_top']+i[1][f'{depth_ref}_bottom'])/2),xycoords='data',
                    xytext=(-180, 0), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->"))
                else: 
                    pass
        
            for j in np.arange(i[1][f'{depth_ref}_top'],i[1][f'{depth_ref}_bottom'],0.2):
                vax.hlines(j,0,15,**perf_kw)

    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            vax.hlines(i[1]['depth'],0,1, **corr_kw)
            if cor_ann:
                try:
                    vax.annotate(f"{i[1]['depth']} - {i[1]['comment']}",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    vax.annotate(f"{i[1]['depth']}",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
     
    
    if 'label' in vsh_kw:
        vax.legend()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        