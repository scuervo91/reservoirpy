import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
from scipy.interpolate import interp1d

def cutoff(log,phi=None,sw=None,vsh=None,elbow=0.9,swaq=0.7,ax=None, plot=True):
    step=np.mean(np.diff(log.index))
    Hc=log.loc[log[sw]<swaq,phi]*(1-log.loc[log[sw]<swaq,sw])*step
    Hct=sum(Hc[~np.isnan(Hc)])
    vshr=np.linspace(0,1,20)
    phir=np.linspace(0,0.4,20)
    swr=np.linspace(0,1,20)
    
    ranges=[phir,swr,vshr]
    symbols=[phi,sw,vsh]
    
    hcl=[]
    co={}
    re=[]
    for (i,s) in enumerate(symbols):
        d=[]
        for j in ranges[i]:
            if i!=0:
                x=np.sum(Hc[log[s]<=j])/Hct
            else:
                x=np.sum(Hc[log[s]>=j])/Hct
            d.append(x)
        hcl.append(d)
        f = interp1d(hcl[i],ranges[i])
        co[s]=f(elbow)
    re.append(co)
    if plot:
        ax = ax or plt.gca()
        ax.plot(phir,hcl[0],color='red',label='Phie CutOff={}'.format(np.round(co[symbols[0]],decimals=2)))
        ax.hlines(elbow,0,1,linestyle='--')
        ax.vlines(co[symbols[0]],0,1,linestyle='--')
        ax.set_xlim([0,1])
        ax.set_ylim([0,1])
        ax.plot(swr,hcl[1],color='blue',label="Sw CutOff={}".format(np.round(co[symbols[1]],decimals=2)))
        ax.vlines(co[symbols[1]],0,1,linestyle='--')
        ax.plot(vshr,hcl[2],color='gold',label='Vsh CutOff={}'.format(np.round(co[symbols[2]],decimals=2)))       
        ax.vlines(co[symbols[2]],0,1,linestyle='--')
        ax.set_xlabel('Phi-Sw-Vsh CutOff') 
        ax.set_ylabel('Norm Hydrocarbon Column')
        ax.set_title('CutOff Estimation HC')
        ax.legend()
        re.append(ax)
    return re

def pickett(rw=0.15,a=1,m=2,n=2,swr=np.linspace(0.1,1,5)):

    Phi_rt1=np.power(rw,1/m)  #Phi at Rt=1. Log(phi)=-mlog(Rt/Rw)
    rts=(a*rw)/(np.power(Phi_rt1,m)*np.power(swr,n))
    Rti=rts*np.power(Phi_rt1,m)/a
    phir=np.logspace(-2,0,num=9)
    RT=[]
    for i in Rti:
        rtr=i*a*np.power(phir,-m)
        RT.append(rtr)
    
    ax=plt.gca()
    for (k,i) in enumerate(RT):
        ax.loglog(i,phir,label='Sw={}'.format(swr[k]))
    ax.set_xticks(np.logspace(-1,3,5))
    ax.set_xticklabels(np.logspace(-1,3,5))
    ax.set_yticks(np.logspace(-2,0,3))
    ax.set_yticklabels(np.logspace(-2,0,3))
    ax.legend()
    ax.set_ylabel('Porosity[]')
    ax.set_xlabel('Resistivity [Ohm m]')
    ax.set_title('Pickett Plot- a={},m={},n={},rw={}'.format(a,m,n,rw))
    return ax

def windland(phirange=None,r35range=None):
    if phirange is None:
        phirange=np.linspace(0.05,0.30,5)
    else:
        phirange=np.asarray(phirange)
    if r35range is None:
        r35range=np.logspace(0,1,5)
    else:
        r35range=np.asarray(r35range)
    
    A=0.735
    B=0.588
    C=0.8641
    K=np.zeros((phirange.shape[0],r35range.shape[0]))
    for (i,v) in enumerate(r35range):
        k=np.power((v*np.power(phirange*100,C)*np.exp(-A)),1/B)
        K[:,i]=k
        
    ax=plt.gca()
    ax.plot(phirange,K)
    ax.set_yscale('log')
    ax.set_xlabel('Phi[]')
    ax.set_ylabel('Permeability [md]')
    ax.legend([i for i in r35range])
    return ax
