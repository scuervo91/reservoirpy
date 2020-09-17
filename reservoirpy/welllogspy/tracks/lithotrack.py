import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np

def lithotrack(df:pd.DataFrame,
               codecols:list,
               percols:list, 
               dtick:bool=False, 
               lims:list=None, 
               codedict: dict=None, 
               fontsize=8, 
               ax=None,
               correlation: pd.DataFrame = None,
               grid_numbers : list = [11,51],
               steps: list  = None,
               corr_kw={}):
    """lithotrack [summary]

    Parameters
    ----------
    df : pd.DataFrame
        [description]
    codecols : list
        [description]
    percols : list
        [description]
    dtick : bool, optional
        [description], by default False
    lims : list, optional
        [description], by default None
    codedict : dict, optional
        [description], by default None
    fontsize : int, optional
        [description], by default 8
    ax : [type], optional
        [description], by default None
    correlation : pd.DataFrame, optional
        [description], by default None
    grid_numbers : list, optional
        [description], by default [11,51]
    steps : list, optional
        [description], by default None
    corr_kw : dict, optional
        [description], by default {}
    """
    lit=ax or plt.gca()
    
    def_corr_kw = {
    'color': 'red',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_corr_kw.items():
        if k not in corr_kw:
            corr_kw[k]=v
    
    df.index.names=['depth']
    df=df.reset_index()
    df=df.loc[(df.index>=lims[0])&(df.index<=lims[1]),:]
    #Create a pivot table concatenating the lithology code names
    mm=pd.DataFrame()
    for (k,v) in enumerate(codecols):    
        m=df.pivot_table(index=['depth'],columns=[v],values=percols[k])
        mm=pd.concat([mm,m],axis=1)

    #Merge in a single dataframe the repeated colnames
    mm=mm.fillna(0)
    lm=pd.DataFrame()
    for i in mm.columns.unique():
        if mm[i].ndim>1:             
            lm[i]=mm[i].max(axis=1)
        elif mm[i].ndim==1:
            lm[i]=mm[i]
    try:
        lm=lm.drop(columns=[0])
    except:
        pass
    lmc=np.cumsum(lm,axis=1)
    
    for i, col in enumerate(lmc.columns):
        lit.fill_betweenx(lmc.index, lmc.iloc[:,i], label=codedict[col], zorder=-i)
        
    if lims==None: #Depth Limits
        lims=[df.index.min(),df.index.max()]

    lit.set_ylim([lims[1],lims[0]])
        
    #Set the vertical grid spacing
    if steps is None:
        mayor_grid = np.linspace(lims[0],lims[1],grid_numbers[0])
        minor_grid = np.linspace(lims[0],lims[1],grid_numbers[1])
    else:
        mayor_grid = np.arange(lims[0],lims[1],steps[0])
        minor_grid = np.arange(lims[0],lims[1],steps[1])
        
    lit.legend()
    lit.set_xlim([0,100])
    lit.set_yticks(mayor_grid)
    lit.set_yticks(minor_grid,minor=True)  
    if dtick==True:
        lit.set_yticklabels(mayor_grid)
    else:
        lit.set_yticklabels([])
    lit.set_xlabel("Lithology")
    lit.xaxis.tick_top()
    lit.xaxis.set_label_position("top")
    lit.tick_params("both",labelsize=fontsize)    

def lithointtrack(depth,lint=None,dtick=False, lims=None, codedict=None, fontsize=8, ax=None):
    lth = ax or plt.gca()
    master=lint[(lint.index>=lims[0])&(lint.index<=lims[1])]
    m=pd.get_dummies(master)
    try: 
        m=m.drop(columns=[0])
    except:
        pass
    for i, col in enumerate(m.columns):
        lth.fill_betweenx(m.index, m.iloc[:,i], label=codedict[col], zorder=-i)
        
    if lims==None: #Depth Limits
        lims=[depth.max(),depth.min()]
        lth.set_ylim(lims)
    else:
        lth.set_ylim([lims[1],lims[0]])
    lth.legend()
    lth.set_xlim([0,1])
    lth.set_yticks(np.linspace(lims[0],lims[1],11))
    lth.set_yticks(np.linspace(lims[0],lims[1],51),minor=True)  
    if dtick==True:
        lth.set_yticklabels(np.linspace(lims[0],lims[1],11))
    else:
        lth.set_yticklabels([])
    lth.set_xlabel("Lithology Interpreted")
    lth.xaxis.tick_top()
    lth.xaxis.set_label_position("top")
    lth.tick_params("both",labelsize=fontsize)