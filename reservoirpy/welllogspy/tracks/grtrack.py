import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import matplotlib as mpl

def grtrack(df: pd.DataFrame,
                gr: (str,list) = None,
                sp: (str,list) = None,
                lims: list = None,
                gr_max: int = 200,
                sp_lim:list = None,
                fm: pd.DataFrame = None,
                units: pd.DataFrame = None,
                perf: pd.DataFrame = None,
                gr_sand_shale: pd.DataFrame = None,
                correlation: pd.DataFrame = None,
                fontsize: int = 8,
                dtick: bool = True,
                sand_flag: pd.Series=None,
                grid_numbers : list = [11,51],
                steps: list  = None,
                legend = False,
                gr_kw={},
                sp_kw={},
                fm_kw={},
                unit_kw={},
                perf_kw={},
                gr_sand_kw={},
                gr_shale_kw={},
                corr_kw={},
                ax=None,
                depth_ref='md',
                gr_colormap: str='autumn',
                sp_colormap: str='gray',
                sp_norm = True
               ):
    """grtrack [summary]

    Parameters
    ----------
    df : pd.DataFrame
        [description]
    gr : [type], optional
        [description], by default None
    sp : [type], optional
        [description], by default None
    lims : list, optional
        [description], by default None
    gr_max : int, optional
        [description], by default 200
    sp_lim : list, optional
        [description], by default None
    fm : pd.DataFrame, optional
        [description], by default None
    units : pd.DataFrame, optional
        [description], by default None
    perf : pd.DataFrame, optional
        [description], by default None
    gr_sand_shale : pd.DataFrame, optional
        [description], by default None
    correlation : pd.DataFrame, optional
        [description], by default None
    fontsize : int, optional
        [description], by default 8
    dtick : bool, optional
        [description], by default True
    sand_flag : pd.Series, optional
        [description], by default None
    grid_numbers : list, optional
        [description], by default [11,51]
    steps : list, optional
        [description], by default None
    legend : bool, optional
        [description], by default False
    gr_kw : dict, optional
        [description], by default {}
    sp_kw : dict, optional
        [description], by default {}
    fm_kw : dict, optional
        [description], by default {}
    unit_kw : dict, optional
        [description], by default {}
    perf_kw : dict, optional
        [description], by default {}
    gr_sand_kw : dict, optional
        [description], by default {}
    gr_shale_kw : dict, optional
        [description], by default {}
    corr_kw : dict, optional
        [description], by default {}
    ax : [type], optional
        [description], by default None
    depth_ref : str, optional
        [description], by default 'md'
    gr_colormap : str, optional
        [description], by default 'autumn'
    sp_colormap : str, optional
        [description], by default 'gray'
    sp_norm : bool, optional
        [description], by default True
    """
    assert isinstance(df,pd.DataFrame)
    assert depth_ref in ['md','tvd','tvdss'], "depth_ref can only be one of ['md','tvd','tvdss']"

    #Create the Axex
    grax= ax or plt.gca()
    
    # Default kwargs for all Lines
    def_gr_kw = {
    'color': 'darkred',
    'linestyle':'-',
    'linewidth': 2
    }
    for (k,v) in def_gr_kw.items():
        if k not in gr_kw:
            gr_kw[k]=v

    def_sp_kw = {
    'color': 'darkblue',
    'linestyle':'-',
    'linewidth': 1
    }    
    for (k,v) in def_sp_kw.items():
        if k not in sp_kw:
            sp_kw[k]=v
    
    def_fm_kw = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 2
    }    
    for (k,v) in def_fm_kw.items():
        if k not in fm_kw:
            fm_kw[k]=v

    def_unit_kw = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 2
    }    
    for (k,v) in def_unit_kw.items():
        if k not in unit_kw:
            unit_kw[k]=v

    def_perf_kw = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 1
    }    
    for (k,v) in def_perf_kw.items():
        if k not in perf_kw:
            perf_kw[k]=v
    
    def_gr_sand_kw = {
    'color': 'gold',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_gr_sand_kw.items():
        if k not in gr_sand_kw:
            gr_sand_kw[k]=v

    def_gr_shale_kw = {
    'color': 'gray',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_gr_shale_kw.items():
        if k not in gr_shale_kw:
            gr_shale_kw[k]=v

    def_corr_kw = {
    'color': 'red',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_corr_kw.items():
        if k not in corr_kw:
            corr_kw[k]=v

 ## Plot the Main Curves      
    #If GammaRay Provided
    depth = df.index if depth_ref=='md' else df[depth_ref]   
    if gr is not None:
        if isinstance(gr,str):
            grax.plot(df[gr],depth,**gr_kw)   #Plotting
        elif isinstance(gr,list):
            cmap = mpl.cm.get_cmap(gr_colormap,len(gr))
            for i,g in enumerate(gr):
                gr_kw['color']=cmap(i)
                grax.plot(df[g],depth,**gr_kw)
        
     #If sp provided   
    if sp is not None:
        spax=grax.twiny()
        if isinstance(sp,str):
            sp_filter = df.loc[df[sp]>-999,sp]
            sp_mean = sp_filter.mean() if sp_norm else 0
            sp_norm = sp_filter - sp_mean
            sp_min = sp_norm.min()
            sp_max = sp_norm.max()
            spax.plot(df[sp]-sp_mean,depth,**sp_kw)   #Plotting
        elif isinstance(sp,list):
            cmap = mpl.cm.get_cmap(sp_colormap,len(sp))
            for i,s in enumerate(sp):
                sp_kw['color']=cmap(i)
                sp_filter = df.loc[df[s]>-999,s]
                sp_mean = sp_filter.mean() if sp_norm else 0
                sp_norm = sp_filter - sp_mean
                sp_min = sp_norm.min()
                sp_max = sp_norm.max()
                spax.plot(df[s]-sp_mean,depth,**sp_kw)

        spax.xaxis.set_label_position("bottom")
        spax.xaxis.tick_bottom()
        spax.set_xlabel('Sp')
       
        if sp_lim is None:
            spax.set_xlim([sp_min,sp_max])
            spax.set_xticks(np.linspace(sp_min,sp_max,4))
 
        else:
            spax.set_xlim(sp_lim)
            spax.set_xticks(np.linspace(sp_lim[0],sp_lim[1],5))
 
        

        
    # Set The lims of depth    
    grax.set_xlim([0,gr_max])           
    if lims==None: #Depth Limits
        lims=[depth.min(),depth.max()]

    grax.set_ylim([lims[1],lims[0]])

    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])

    #Set Gridding and ticks
    grax.set_xlabel('GammaRay')
    grax.tick_params("both",labelsize=fontsize)
    grax.grid(True,linewidth=1.0)
    grax.grid(True,which='minor', linewidth=0.5)
    plt.setp(grax.get_yticklabels(),visible=True)
    grax.set_yticks(minor_grid,minor=True)
    grax.set_yticks(mayor_grid)
    if dtick==True:
        grax.set_yticklabels(mayor_grid)
    else:
        grax.set_yticklabels([])
    grax.set_xticks(np.linspace(0,gr_max,4))
    grax.xaxis.set_label_position("top")
    grax.xaxis.tick_top()
    
    #Add Formations tops
    if fm is not None:
        fm_ann = fm_kw.pop('ann',False)
        fm_ann_fontsize = fm_kw.pop('fontsize',8)
        for i in fm.iterrows():
            if depth_ref == 'tvdss':
                if i[1][f'{depth_ref}_top'] >= lims[0] or i[1][f'{depth_ref}_top'] <= lims[1]:
                    continue
            else:
                if i[1][f'{depth_ref}_top'] <= lims[0] or i[1][f'{depth_ref}_top'] >= lims[1]:
                    continue
            grax.hlines([i[1][f'{depth_ref}_top']],0,gr_max, **fm_kw)
            if fm_ann:
                grax.annotate(f"Top of {i[0]}",xy=(gr_max-3,i[1][f'{depth_ref}_top']-2),
                                xycoords='data',horizontalalignment='right', fontsize=fm_ann_fontsize,
                                bbox={'boxstyle':'round', 'fc':'0.8'})

    #Add units tops
    if units is not None:
        unit_ann = unit_kw.pop('ann',False)
        unit_ann_fontsize = unit_kw.pop('fontsize',8)
        for i in units.iterrows():
            if depth_ref == 'tvdss':
                if i[1][f'{depth_ref}_top'] >= lims[0] or i[1][f'{depth_ref}_top'] <= lims[1]:
                    continue
            else:
                if i[1][f'{depth_ref}_top'] <= lims[0] or i[1][f'{depth_ref}_top'] >= lims[1]:
                    continue
            grax.hlines([i[1][f'{depth_ref}_top']],0,gr_max, **unit_kw)
            if unit_ann:
                grax.annotate(f"Top of {i[0]}",xy=(gr_max-3,i[1][f'{depth_ref}_top']-2),
                                xycoords='data',horizontalalignment='right', fontsize=unit_ann_fontsize,
                                bbox={'boxstyle':'round', 'fc':'0.8'})
                          

     #Add Interval Perforated
    if perf is not None:
        perf_ann = perf_kw.pop('ann',False)
        perf_ann_fontsize = perf_kw.pop('fontsize',8)
        perf_ann_xytext = perf_kw.pop('xytext',(-180,0))
        for i in perf.iterrows():
            if depth_ref == 'tvdss':
                if i[1][f'{depth_ref}_top'] >= lims[0] or i[1][f'{depth_ref}_top'] <= lims[1]:
                    continue
            else:
                if i[1][f'{depth_ref}_top'] <= lims[0] or i[1][f'{depth_ref}_top'] >= lims[1]:
                    continue
            if perf_ann:
                try:
                    grax.annotate(f"Top:{i[1][f'{depth_ref}_top']} \nBottom:{i[1][f'{depth_ref}_bottom']} \nNote:{i[1]['comment']}",
                                  xy=(0,(i[1][f'{depth_ref}_top']+i[1][f'{depth_ref}_bottom'])/2),xycoords='data',
                                  xytext=(-180, 0), textcoords='offset points',
                                  arrowprops=dict(arrowstyle="->"))
                except:
                    grax.annotate(f"Top:{i[1][f'{depth_ref}_top']} \nBottom:{i[1][f'{depth_ref}_bottom']}",
                    xy=(0,(i[1][f'{depth_ref}_top']+i[1][f'{depth_ref}_bottom'])/2),xycoords='data',
                    xytext=perf_ann_xytext, textcoords='offset points', fontsize = perf_ann_fontsize,
                    arrowprops=dict(arrowstyle="->"))
                else: 
                    pass
            if depth_ref == 'tvdss':
                for j in np.arange(i[1][f'{depth_ref}_bottom'],i[1][f'{depth_ref}_top'],0.5):
                    grax.hlines(j,0,15,**perf_kw)
            else:
                for j in np.arange(i[1][f'{depth_ref}_top'],i[1][f'{depth_ref}_bottom'],0.5):
                    grax.hlines(j,0,15,**perf_kw)

    #Add Sand Gamma Ray line      
    if gr_sand_shale is not None:
        for i in gr_sand_shale.iterrows():
            try:
                grax.vlines(i[1]['gr_sand'],i[1][f'{depth_ref}_top'],i[1][f'{depth_ref}_bottom'],**gr_sand_kw)
            except:
                pass
            
            try:
                grax.vlines(i[1]['gr_shale'],i[1][f'{depth_ref}_top'],i[1][f'{depth_ref}_bottom'],**gr_shale_kw)
            except:
                pass
             
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            if depth_ref == 'tvdss':
                if i[1]['depth'] >= lims[0] or i[1]['depth'] <= lims[1]:
                    continue
            else:
                if i[1]['depth'] < lims[0] or i[1]['depth'] > lims[1]:
                    continue
            grax.hlines(i[1]['depth'],0,gr_max, **corr_kw)
            if cor_ann:
                try:
                    grax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(gr_max-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    grax.annotate(f"{i[1]['depth']}",xy=(gr_max-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                
    
    if legend:
        grax.legend() 