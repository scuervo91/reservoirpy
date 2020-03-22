import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np


def permtrack(depth, k=None, lims=None,
            dtick=False,
            krange=[0.2,20000],ax=None,
            fontsize=8,**kw):

    kax=ax or plt.gca()
    defkwa = {
    'linestyle':'-',
    'linewidth': 1
    }
    
    for (ky,v) in defkwa.items():
        if ky not in kw:
            kw[ky]=v
    if k is not None:
        kax.plot(k,depth,**kw)
    
    if lims==None: #Depth Limits
        lims=[depth.max(),depth.min()]
        kax.set_ylim(lims)
    else:
        kax.set_ylim([lims[1],lims[0]])
        
    kax.set_xscale("log")
    kax.set_xlim(krange)
    ticks=np.round(np.power(10,np.linspace(np.log10(krange[0]),np.log10(krange[1]),np.log10(krange[1]/krange[0])+1)),decimals=1)
    kax.set_xticks(ticks)
    kax.set_xticklabels(ticks)
    kax.set_xlabel("Permeability [md]")
    kax.xaxis.tick_top()
    kax.xaxis.set_label_position("top")
    kax.tick_params("both",labelsize=fontsize)
    kax.xformatter = 'auto'
    kax.set_yticks(np.linspace(lims[0],lims[1],11))
    kax.set_yticks(np.linspace(lims[0],lims[1],51),minor=True)
    if dtick==True:
        kax.set_yticklabels(np.linspace(lims[0],lims[1],11))
    else:
        kax.set_yticklabels([])
    kax.grid(True,linewidth=1.0)
    kax.grid(True,which='minor', linewidth=0.5)
    if 'label' in kw:
        kax.legend()