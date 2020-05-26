import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import numpy as np

def dntrack(df: pd.DataFrame,
            rho: (list,str) = None,
            ntr: (list,str) = None,
            lims: list = None,
            lime: bool = False,
            dtick: bool =False,
            fill: bool =True,
            fontsize: int=8,
            grid_numbers : list = [11,51],
            steps: list  = None,
            correlation: pd.DataFrame = None,
            rho_kw:dict={},
            ntr_kw:dict={},
            corr_kw:dict={},
            ax=None,
            rho_colormap:str='hot',
            ntr_colormap:str='winter',
            depth_ref:str='md'
           ):
    """
    Parameters
    ----------
    df : pd.DataFrame
        Series with Depth Array
    rho : str, optional
        Series with Density Array. The default is None.
    ntr : str, optional
        Series with Neutron Array. The default is None.
    lims : list, optional
        list that contains the depth limits to plot. e.g lims = [5800,6500]. The default is None.
    lime : bool, optional
        Define the logs calibration depending if limestone or sandstone. The default is False.
    dtick : bool, optional
        Bool If show the depth ticks.. The default is False.
    fill : bool, optional
        Bool if fill the logs cross. The default is True.
    fontsize : int, optional
       Ticks Fontsize. The default is 8.
    grid_numbers : list, optional
        List indicates the number of the grids to show in the depth 
        axis for mayor and minor. The default is [11,51].
    steps : list, optional
        List indicates the step size of the grids to show in the depth 
        axis for mayor and minor. If 'steps' is provided 'grid_numbers'
        is not used. The default is None.
    correlation : plot correlate lines at specific depth. he default is None
    rho_kw : TYPE, optional
        kwargs for Density curve. The default is {}.
    ntr_kw : TYPE, optional
        kwargs for Neutron curve.. The default is {}.
    corr_kw : TYPE, optional
        kwargs for Neutron curve.. The default is {}.
    ax : TYPE, optional
        Ax to plot. The default is None.

    Returns
    -------
    Figure

    """
    assert isinstance(df,pd.DataFrame)
    assert depth_ref in ['md','tvd','tvdss'], "depth_ref can only be one of ['md','tvd','tvdss']"

    #Set Axes
    dax=ax or plt.gca()
    nax=dax.twiny()

    # Default kwargs for rho and ntr lines
    def_rho_kw = {
    'color': 'darkred',
    'linestyle':'-',
    'linewidth': 2
    }
    for (k,v) in def_rho_kw.items():
        if k not in rho_kw:
            rho_kw[k]=v

    def_ntr_kw = {
    'color': 'darkblue',
    'linestyle':'-',
    'linewidth': 1
    }    
    for (k,v) in def_ntr_kw.items():
        if k not in ntr_kw:
            ntr_kw[k]=v
            
    def_corr_kw = {
    'color': 'red',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_corr_kw.items():
        if k not in corr_kw:
            corr_kw[k]=v
  
    #Set type of sync between Neutron GammaRay
    if lime==True:
        d=2.71
    else:
        d=2.65

    m=(d-1.9)/(0-0.45)
    b=-m*0.45+1.9
    rholim=-0.15*m+b

    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
    
    depth = df.index if depth_ref=='md' else df[depth_ref]
    #Set Density Axes
    if rho is not None:
        if isinstance(rho,str):
            dax.plot(df[rho],depth,**rho_kw)   #Plotting
        elif isinstance(rho,list):
            cmap = mpl.cm.get_cmap(rho_colormap,len(rho))
            for i,r in enumerate(rho):
                rho_kw['color']=cmap(i)
                dax.plot(df[r],depth,**rho_kw)

        #Set the gridding and ticks
        dax.set_xlabel("Density [g/cc]")
        dax.set_xticks(np.linspace(1.9,rholim,4))
        dax.set_xlim([1.9,rholim])
        dax.tick_params("both",labelsize=fontsize)
        dax.grid(True,linewidth=1.0)
        dax.grid(True,which='minor', linewidth=0.5)
        dax.set_yticks(minor_grid,minor=True)
        dax.set_yticks(mayor_grid)
        if dtick==True:
            dax.set_yticklabels(mayor_grid)
        else:
            dax.set_yticklabels([])

    
    #Set neutron axes
    if ntr is not None:
        if isinstance(ntr,str):
            nax.plot(df[ntr],depth,**ntr_kw)   #Plotting
        elif isinstance(ntr,list):
            cmap = mpl.cm.get_cmap(ntr_colormap,len(ntr))
            for i,r in enumerate(ntr):
                ntr_kw['color']=cmap(i)
                nax.plot(df[r],depth,**ntr_kw)
    
        nax.set_xlabel("Neutron [v/v]")
        nax.set_xticks(np.linspace(0.45,-0.15,4))
        nax.set_xlim([0.45,-0.15])
        nax.tick_params("both",labelsize=fontsize)
        nax.set_yticks(minor_grid,minor=True)
        nax.set_yticks(mayor_grid)
        if dtick==True:
            nax.set_yticklabels(mayor_grid)
        else:
            nax.set_yticklabels([])


    if lims==None: #Depth Limits
        lims=[depth.min(),depth.max()]
        
    dax.set_ylim([lims[1],lims[0]])
    #Convert the Neutron values to Density Values in order to fill the cross Density-Neutron
    #When the track is callibrated for sandstone use m=-1.666667 and b=2.65
    #When the track is callibrated for limestone use m=-1.8 and b=2.71
    if (ntr is not None) & (ntr is not None):
        NtrTorho=df[ntr]*m+b
        ntrrho=NtrTorho.values.ravel()

        if fill==True:
            dax.fill_betweenx(depth,df[rho],ntrrho,where=(df[rho] < ntrrho),color="red")   
            
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            dax.hlines(i[1]['depth'],0,rholim, **corr_kw)
            if cor_ann:
                try:
                    dax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(rholim-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    dax.annotate(f"{i[1]['depth']}",xy=(rholim-3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})