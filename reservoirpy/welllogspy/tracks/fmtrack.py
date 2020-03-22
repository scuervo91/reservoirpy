import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np

def fmtrack(fm=None,lims=None, top='top',bottom='bottom',name='name',ax=None):
    
    fax=ax or plt.gca()
        
    fm['mid']=fm.mean(axis=1)
    if lims==None: #Depth Limits
        lims=[fmf.loc[:,top].max(),fmf.loc[:,top].max()]
    fax.set_ylim([lims[1],lims[0]])
    fmf=fm[(fm['mid']<=lims[1]) & (fm['mid']>=lims[0])]
    m=fmf.shape[0]
    
    for i in fmf.iterrows():
        fax.axhspan(i[1][top],i[1][bottom],xmin=0,xmax=1,alpha=0.2)
        fax.hlines([i[1][top],i[1][bottom]],0,1, linewidth=1.0)
        
    fax.set_yticks(fmf.mid)
    fax.set_yticklabels(fmf[name])