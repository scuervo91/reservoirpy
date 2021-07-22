import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import matplotlib as mpl
from typing import Union, Sequence

def grtrack(df: pd.DataFrame,
                gr: Union[str,list] = None,
                sp: Union[str,list] = None,
                lims: Sequence[float] = None,
                gr_max: int = 200,
                sp_lim:list = None,
                formations: pd.DataFrame = None,
                units: pd.DataFrame = None,
                perf: pd.DataFrame = None,
                gr_sand_shale: pd.DataFrame = None,
                correlation: pd.DataFrame = None,
                fontsize: int = 8,
                dtick: bool = True,
                grid_numbers : Sequence[float] = [11,51],
                steps: Sequence[float]  = None,
                gr_steps=4,
                legend:bool = False,
                gr_kw:dict={},
                sp_kw:dict={},
                formation_kw:dict={},
                unit_kw:dict={},
                perf_kw:dict={},
                gr_sand_kw:dict={},
                gr_shale_kw:dict={},
                corr_kw:dict={},
                ax=None,
                depth_ref:str='md',
                gr_colormap: str='autumn',
                sp_colormap: str='gray',
                sp_norm:bool = False,
                sp_baseline:float = 0,
                sp_fill:bool = False
               ):
    """grtrack [Create a Gamma Ray & SP track with matplotlib as a backend]

    Parameters
    ----------
    df : pd.DataFrame
        [DataFrame with the data to plot indexed by the depth]
    gr : Union[str,list], optional
        [Gamma Ray Columns name in the df DataFrame], by default None
    sp : Union[str,list], optional
        [Gamma Ray Columns name in the df DataFrame], by default None
    lims : Sequence[float], optional
        [Interval depth to plot. List of numbers [Top, Bottom]], by default None
    gr_max : int, optional
        [Set max Gr interval to plot], by default 200
    sp_lim : list, optional
        [Set sp lim], by default None
    fm : pd.DataFrame, optional
        [Add Tops depth], by default None
    units : pd.DataFrame, optional
        [Add Units depth], by default None
    perf : pd.DataFrame, optional
        [Add perforations depth], by default None
    gr_sand_shale : pd.DataFrame, optional
        [Add Gr sand a shale depth to plot], by default None
    correlation : pd.DataFrame, optional
        [Add correlations], by default None
    fontsize : int, optional
        [Plot Fontsize], by default 8
    dtick : bool, optional
        [If True show the depth tick. ], by default True
    grid_numbers : Sequence[float], optional
        [Number of mayor and minor horizontal grids], by default [11,51]
    steps : Sequence[float], optional
        [Set size of horizontal grids], by default None
    legend : bool, optional
        [If true show the curves legend], by default False
    gr_kw : dict, optional
        [Add key arguments that modify the gamma ray curve], by default {}
    sp_kw : dict, optional
        [Add key arguments that modify the sp curve], by default {}
    formation_kw : dict, optional
        [Add key arguments that modify the formation curve], by default {}
    unit_kw : dict, optional
        [Add key arguments that modify the units curve], by default {}
    perf_kw : dict, optional
        [Add key arguments that modify the perforations curve], by default {}
    gr_sand_kw : dict, optional
        [Add key arguments that modify the Sand vertical curve], by default {}
    gr_shale_kw : dict, optional
        [Add key arguments that modify the Shale curve], by default {}
    corr_kw : dict, optional
        [Add key arguments that modify the correlation curven], by default {}
    ax : [type], optional
        [Matplotlib axes], by default None
    depth_ref : str, optional
        [Set the depth reference: md, tvd, tvdss], by default 'md'
    gr_colormap : str, optional
        [Set colormap when there are multiple gammaray curves], by default 'autumn'
    sp_colormap : str, optional
        [Set colormap when there are multiple sp curves], by default 'gray'
    sp_norm : bool, optional
        [Normalize sp curve], by default False
    sp_baseline : float, optional
        [Plot the sp baseline], by default 0
    sp_fill : bool, optional
        [Fill the sp curve], by default False
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
    'linewidth': 1,
    'sp_left_color':(0.9,0.57,0.2,0.2),
    'sp_right_color':(0.28,0.9,0.79,0.2)
    }    
    for (k,v) in def_sp_kw.items():
        if k not in sp_kw:
            sp_kw[k]=v
    
    def_formation_kw = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 2,
    'fill':False,
    'cmap':'jet',
    'alpha':0.15
    }    
    for (k,v) in def_formation_kw.items():
        if k not in formation_kw:
            formation_kw[k]=v

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
            sp_left_color = sp_kw.pop('sp_left_color',(0.9,0.57,0.2,0.2))
            sp_right_color = sp_kw.pop('sp_right_color',(0.28,0.9,0.79,0.2))
            sp_filter = df.loc[df[sp]>-999,sp]
            sp_mean = sp_filter.mean() if sp_norm else sp_baseline
            sp_norm = sp_filter - sp_mean
            sp_min = sp_norm.min()
            sp_max = sp_norm.max()
            spax.plot(df[sp]-sp_mean,depth,**sp_kw)   #Plotting
            if sp_fill==True:
                spax.fill_betweenx(depth,df[sp]-sp_mean,0,where=(df[sp]-sp_mean < 0),color=sp_left_color)
                spax.fill_betweenx(depth,df[sp]-sp_mean,0,where=(df[sp]-sp_mean > 0),color=sp_right_color)

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
        if depth_ref=='tvdss':
            mayor_grid = np.arange(lims[1],lims[0],steps[0])
            minor_grid = np.arange(lims[1],lims[0],steps[1])
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
    grax.set_xticks(np.linspace(0,gr_max,gr_steps))
    grax.xaxis.set_label_position("top")
    grax.xaxis.tick_top()
    
    #Add formations tops
    if formations is not None:
        formation_ann = formation_kw.pop('ann',False)
        formation_ann_fontsize = formation_kw.pop('fontsize',8)
        formation_fill = formation_kw.pop('fill',False)
        formation_cmap = formation_kw.pop('cmap','jet')
        formation_alpha = formation_kw.pop('alpha',0.2)

        if depth_ref == 'tvdss':
            formations = formations.loc[(formations[f'{depth_ref}_top']<= lims[0])&(formations[f'{depth_ref}_top']>= lims[1]),:]
        else:
            formations = formations.loc[(formations[f'{depth_ref}_top']>= lims[0])&(formations[f'{depth_ref}_top']<= lims[1]),:]
        formation_cmap_shape = formation_kw.pop('cmap_shape',formations.shape[0])
        
        if 'color_id' not in formations.columns:
            formations['color_id'] = np.arange(formations.shape[0])
        for i,r in formations.iterrows():
            # Fill with color between the top and bottom of each formation
            if formation_fill:
                if depth_ref == 'tvdss':
                    fill_top = lims[0] if r[f'{depth_ref}_top'] >= lims[0] else r[f'{depth_ref}_top']
                    fill_bottom = lims[1] if r[f'{depth_ref}_bottom'] <= lims[1] else r[f'{depth_ref}_bottom']
                else:
                    fill_top = lims[0] if r[f'{depth_ref}_top'] <= lims[0] else r[f'{depth_ref}_top']
                    fill_bottom = lims[1] if r[f'{depth_ref}_bottom'] >= lims[1] else r[f'{depth_ref}_bottom']
                
                grax.fill_between([0,gr_max],fill_top,fill_bottom,color=mpl.cm.get_cmap(formation_cmap,formation_cmap_shape)(r['color_id'])[:3]+(formation_alpha,))

            grax.hlines([r[f'{depth_ref}_top']],0,gr_max, **formation_kw)
            if formation_ann:
                grax.annotate(f"Top of {i}",xy=(gr_max-3,r[f'{depth_ref}_top']-2),
                                xycoords='data',horizontalalignment='right', fontsize=formation_ann_fontsize,
                                bbox={'boxstyle':'round', 'fc':'0.8'})

    #Add units tops
    if units is not None:
        unit_ann = unit_kw.pop('ann',False)
        unit_ann_fontsize = unit_kw.pop('fontsize',8)
        unit_fill = unit_kw.pop('fill',False)
        unit_cmap = unit_kw.pop('cmap','jet')
        unit_alpha = unit_kw.pop('alpha',0.2)

        if depth_ref == 'tvdss':
            units = units.loc[(units[f'{depth_ref}_top']<= lims[0])&(units[f'{depth_ref}_top']>= lims[1]),:]
        else:
            units = units.loc[(units[f'{depth_ref}_top']>= lims[0])&(units[f'{depth_ref}_top']<= lims[1]),:]
        unit_cmap_shape = unit_kw.pop('cmap_shape',units.shape[0])
        
        if 'color_id' not in units.columns:
            units['color_id'] = np.arange(units.shape[0])
        for i,r in units.iterrows():
            # Fill with color between the top and bottom of each formation
            if unit_fill:
                if depth_ref == 'tvdss':
                    fill_top = lims[0] if r[f'{depth_ref}_top'] >= lims[0] else r[f'{depth_ref}_top']
                    fill_bottom = lims[1] if r[f'{depth_ref}_bottom'] <= lims[1] else r[f'{depth_ref}_bottom']
                else:
                    fill_top = lims[0] if r[f'{depth_ref}_top'] <= lims[0] else r[f'{depth_ref}_top']
                    fill_bottom = lims[1] if r[f'{depth_ref}_bottom'] >= lims[1] else r[f'{depth_ref}_bottom']
                
                grax.fill_between([0,gr_max],fill_top,fill_bottom,color=mpl.cm.get_cmap(unit_cmap,unit_cmap_shape)(r['color_id'])[:3]+(unit_alpha,))

            grax.hlines([r[f'{depth_ref}_top']],0,gr_max, **unit_kw)
            if unit_ann:
                grax.annotate(f"Top of {i}",xy=(gr_max-3,r[f'{depth_ref}_top']-2),
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
        cor_ann_fontsize = corr_kw.pop('fontsize',8)
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
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'},
                                 fontsize = cor_ann_fontsize)
                except:
                    grax.annotate(f"{i[1]['depth']}",xy=(gr_max-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'},
                                 fontsize = cor_ann_fontsize)
                
    
    if legend:
        grax.legend() 