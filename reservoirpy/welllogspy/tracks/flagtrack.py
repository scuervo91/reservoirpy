import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np

def flagtrack(df: pd.DataFrame,
              sand: str = None,
              res: str = None,
              pay: str = None,
              ax=None,
              lims: list = None,
              fontsize:int =8,
              correlation: pd.DataFrame = None,
              legend:bool = True,
              corr_kw={},
              sand_kw={},
              res_kw={},
              pay_kw={}
              ):
    
    #Create the axes
    fax=ax or plt.gca()
    
    #Default kwargs for lines
    def_sand_kw = {
    'linestyle':'-',
    'linewidth': 0.1
    }
    
    for (k,v) in def_sand_kw.items():
        if k not in sand_kw:
            sand_kw[k]=v
            
    def_res_kw = {
    'linestyle':'-',
    'linewidth': 0.1
    }
    
    for (k,v) in def_res_kw.items():
        if k not in res_kw:
            res_kw[k]=v
            
    def_pay_kw = {
    'linestyle':'-',
    'linewidth': 0.1
    }
    
    for (k,v) in def_pay_kw.items():
        if k not in pay_kw:
            pay_kw[k]=v
    
    def_corr_kw = {
    'color': 'red',
    'linestyle':'--',
    'linewidth': 2
    }    
    for (k,v) in def_corr_kw.items():
        if k not in corr_kw:
            corr_kw[k]=v
            

    
    # Plot main Lines
    if sand is not None:
        fax.plot(df[sand]*0.33, df.index,**sand_kw)
        fax.fill_betweenx(df.index,0,df[sand]*0.33, color='yellow',label='Sand')
        
    if res is not None:
        fax.plot(df[res]*0.66, df.index,**res_kw)
        fax.fill_betweenx(df.index,0.33,df[res]*0.66, where=(df[res]*0.66>0.33), color='orange',label='Res')
        
    if pay is not None:
        fax.plot(df[pay], df.index,**pay_kw)
        fax.fill_betweenx(df.index,0.66,df[pay], where=(df[pay]>0.66), color='green',label='pay')

    # Set The lims of depth    
    fax.set_xlim([0,1])           
    if lims==None: #Depth Limits
        lims=[df.index.min(),df.index.max()]

    fax.set_ylim([lims[1],lims[0]])
        
    #Set the vertical grid spacing
        
    #Set Gridding and ticks
    fax.set_xlabel('Flags')
    fax.tick_params("both",labelsize=fontsize)
    fax.xaxis.set_label_position("top")
    fax.set_yticklabels([])
    fax.set_xticklabels([])    
    fax.set_yticks([],minor=True)
    fax.set_yticks([])
    fax.set_xticks([],minor=True)
    fax.set_xticks([])
    #fax.axis('off')
    fax.spines["top"].set_visible(False)
    fax.spines["right"].set_visible(False)
    fax.spines["bottom"].set_visible(False)
    fax.spines["left"].set_visible(False)
    
    #Add Correlation Line
    if correlation is not None:
        cor_ann = corr_kw.pop('ann',False)
        for i in correlation.iterrows():
            fax.hlines(i[1]['depth'],0,1, **corr_kw)
            if cor_ann:
                try:
                    fax.annotate(f"{i[1]['depth']} - {i[1]['comment']} ",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
                except:
                    fax.annotate(f"{i[1]['depth']}",xy=(1-0.3,i[1]['depth']-1),
                                 xycoords='data',horizontalalignment='right',
                                 bbox={'boxstyle':'roundtooth', 'fc':'0.8'})
    
    

    if legend:
        fax.legend()