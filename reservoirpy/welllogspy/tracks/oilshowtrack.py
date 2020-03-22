import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np

def oilshowtrack(depth, oilshow=None, lims=None,
            dtick=False,
            fill=True, ax=None,
            fontsize=8,**kw):

    oax=ax or plt.gca()
    
    defkwa = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 1
    }
    
    for (k,v) in defkwa.items():
        if k not in kw:
            kw[k]=v
    
    if oilshow is not None:
        oax.plot(oilshow,depth,**kw)
    
    if lims==None: #Depth Limits
        lims=[depth.max(),depth.min()]
        oax.set_ylim(lims)
    else:
        oax.set_ylim([lims[1],lims[0]])
    
    oax.set_xlim([0,1])
    oax.set_xlabel("OilShow")
    oax.set_xticks(np.linspace(0,1,4))
    oax.set_xticklabels(np.round(np.linspace(0,1,4),decimals=2))
    oax.xaxis.tick_top()
    oax.xaxis.set_label_position("top")
    oax.tick_params("both",labelsize=fontsize)
    oax.set_yticks(np.linspace(lims[0],lims[1],11))
    oax.set_yticks(np.linspace(lims[0],lims[1],51),minor=True)        
    oax.grid(True,linewidth=1.0)
    oax.grid(True,which='minor', linewidth=0.5)
    if dtick==True:
        oax.set_yticklabels(np.linspace(lims[0],lims[1],11))
    else:
        oax.set_yticklabels([])
    if fill==True:
        oax.fill_betweenx(depth,0,oilshow,color="green")
    if 'label' in kw:
        oax.legend()            