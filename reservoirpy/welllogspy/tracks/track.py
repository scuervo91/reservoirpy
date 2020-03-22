import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np

def track(depth, curve=None, lims=None, clims=None,
            dtick=False,
            fill=True,
            scale='linear',
            curvetitle='Track', ax=None,
            fontsize=8,**kw):

    tax=ax or plt.gca()
    
    defkwa = {
    'color': 'black',
    'linestyle':'-',
    'linewidth': 1
    }
    
    for (k,v) in defkwa.items():
        if k not in kw:
            kw[k]=v
    
    if curve is not None:
        tax.plot(curve,depth,**kw)
    
    if lims==None: #Depth Limits
        lims=[depth.max(),depth.min()]
        tax.set_ylim(lims)
    else:
        tax.set_ylim([lims[1],lims[0]])


    tax.set_xlabel(curvetitle)
    tax.set_xscale(scale)
    if clims==None: 
        clims=[curve.min()*0.9,curve.max()*1.1]

    tax.set_xlim(clims)

    tax.set_xticks(np.linspace(clims[0],clims[1],4))
    tax.set_xticklabels(np.round(np.linspace(clims[0],clims[1],4),decimals=2))
    tax.xaxis.tick_top()
    tax.xaxis.set_label_position("top")
    tax.tick_params("both",labelsize=fontsize)
    tax.set_yticks(np.linspace(lims[0],lims[1],11))
    tax.set_yticks(np.linspace(lims[0],lims[1],51),minor=True)        
    tax.grid(True,linewidth=1.0)
    tax.grid(True,which='minor', linewidth=0.5)
    if dtick==True:
        tax.set_yticklabels(np.linspace(lims[0],lims[1],11))
    else:
        tax.set_yticklabels([])
    if 'label' in kw:
        tax.legend() 