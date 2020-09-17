import numpy as np
import pandas as pd
import scipy.stats as st


def resources(log, phi=None, sw=None, cutoff=None,a=1,h=1,beta=1.1, fluid='oil',m=1000,seed=1706):
    """resources [summary]

    Parameters
    ----------
    log : [type]
        [description]
    phi : [type], optional
        [description], by default None
    sw : [type], optional
        [description], by default None
    cutoff : [type], optional
        [description], by default None
    a : int, optional
        [description], by default 1
    h : int, optional
        [description], by default 1
    beta : float, optional
        [description], by default 1.1
    fluid : str, optional
        [description], by default 'oil'
    m : int, optional
        [description], by default 1000
    seed : int, optional
        [description], by default 1706

    Returns
    -------
    [type]
        [description]
    """
    np.random.seed(seed=seed)
    # Define coefficients depending on the fluid
    if fluid == 'oil':
        c=7.758  #Mbbl
        name='OOIP'
    elif fluid == 'gas':
        c=0.000043560 #Bscf
        name='OGIP'
    r=pd.DataFrame()
    for p in phi:
        phieh = np.histogram(log.loc[log[p]>cutoff[0],p])
        phie_dist = st.rv_histogram(phieh)
        phie_random=phie_dist.rvs(size=m)
        for s in sw:
            swh = np.histogram(log.loc[log[s]<cutoff[1],s])
            sw_dist = st.rv_histogram(swh)
            sw_random=sw_dist.rvs(size=m)
            #Original hydrocarbon in place
            orip = c*a*h*phie_random * (1-sw_random)*(1/beta)
            orip = pd.DataFrame({name:orip})
            orip['PhieCurve']=p
            orip['SwCurve']=s
            orip['Phie']=phie_random
            orip['Sw']=sw_random
            r = r.append(orip)
    return r.reset_index(drop=True)