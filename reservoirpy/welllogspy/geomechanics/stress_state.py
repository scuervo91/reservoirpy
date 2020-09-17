import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def coulomb(mu):
    """coulomb [summary]

    Parameters
    ----------
    mu : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    c = np.power(np.sqrt(np.power(mu,2) + 1) + mu, 2)
    return c

def shmin(sv,pp,mu):
    """shmin [summary]

    Parameters
    ----------
    sv : [type]
        [description]
    pp : [type]
        [description]
    mu : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    smin = ((sv - pp) / coulomb(mu)) + pp
    return smin

def shmax(sv,pp,mu):
    """shmax [summary]

    Parameters
    ----------
    sv : [type]
        [description]
    pp : [type]
        [description]
    mu : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    smax = coulomb(mu) * (sv - pp) + pp
    return smax 

def zobackogram(sv,pp,mu, ax=None):
    """zobackogram [summary]

    Parameters
    ----------
    sv : [type]
        [description]
    pp : [type]
        [description]
    mu : [type]
        [description]
    ax : [type], optional
        [description], by default None
    """
    ax = ax or plt.gca()  

    smin = shmin(sv,pp,mu)
    smax = shmax(sv,pp,mu)

    lines = []
    #45° line. slope=1
    L1 = pd.DataFrame({'shmin':[0,smax*1.1],'shmax':[0,smax*1.1]})
    lines.append(L1)

    #Shmin Line
    L2 = pd.DataFrame({'shmin':[smin,smin],'shmax':[smin,sv]})
    lines.append(L2)

    #Shmax Horizonal line
    L3 = pd.DataFrame({'shmin':[sv,smax],'shmax':[smax,smax]})
    lines.append(L3)

    #Shmax diagonal Line
    L4 = pd.DataFrame({'shmin':[smin,sv],'shmax':[sv,smax]})
    lines.append(L4)
    #Internal Lines
    L5 = pd.DataFrame({'shmin':[smin,sv,sv],'shmax':[sv,sv,smax]})
    lines.append(L5)

    co = ['black','darkblue','darkred','darkgreen','grey']
    ls = ['-','-','-','-','--',]
    wd = [3,2,2,2,1]
    lb = ['45° Line','Shmin Line','Shmax line','Shmax Strike slip','limits']

    for i, v in enumerate(lines):
        ax.plot(v['shmin'],v['shmax'],color = co[i], linestyle = ls[i], linewidth = wd[i])
    ax.set_xlabel('Shmin')
    ax.set_ylabel('Shmax')
    ax.legend(lb)
    ax.set_title('Stress State')