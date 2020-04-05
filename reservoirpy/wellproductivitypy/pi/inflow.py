import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



def oil_j(mu,bo,h=None,k=None,kh=None,re=1490,rw=0.58,s=0):
    if kh is None:
        kh = k*h
    j=(0.00708*kh)/(mu*bo*(np.log(re/rw)-0.75+s))
    return j


def oil_inflow_curve(pr,j,pb,n=10):
    """
    Estimate the Pressure and Flow curves given a Productivity Index parameters

    Attributes:
        pr:     Reservoir Pressure [psi]
        j:      Productivity index [bbl/d*psi]
        pb:     Bubble Pressure [psi]
        n:      Number of steps in the return DataFrame
    
    Return:
        pi:     DataFrame with two columns Pressure [p] and Flow [q]-> Pandas DataFrame
        aof:    Absolute open flow -> Number
    """
    p = np.linspace(0,pr,n)
    q = np.zeros(n)

    if pr <= pb:
        q = (pr * j / 1.8) * (1 - 0.2 *(p / pr) - 0.8 * np.power(p / pr,2))
    else:
        q[p > pb] = j * (pr - p[p > pb])
        q[p <= pb] = (j * (pr - pb)) + (pb * j / 1.8) * (1 - 0.2 *(p[p <= pb] / pb) - 0.8 * np.power(p[p <= pb] / pb,2))
    
    data = pd.DataFrame({'p':p,'q':q})
    aof = data['q'].max()

    return data, aof

class oil_inflow:

    def __init__(self,pr,j,pb,n=10):
        self.pr = pr
        self.j = j
        self.pb = pb 
        self.curve,self.aof = oil_inflow_curve(pr,j,pb,n=n)
        self._pwf_to_flow = interp1d(self.curve['p'],self.curve['q'])
        self._flow_to_pwf = interp1d(self.curve['q'],self.curve['p'])

    def __str__(self):
        return f"""Oil Inflow: 
            Reservoir Pressure: {self.pr} psi 
            Productivity Index: {self.j} bbl/d*psi 
            Bubble Point Pressure: {self.pb} psi  
            AOF = {self.aof} bbl/d 
                """
  
    def __repr__(self):
        return f"""Oil Inflow: 
            Reservoir Pressure: {self.pr} psi 
            Productivity Index: {self.j} bbl/d*psi 
            Bubble Point Pressure: {self.pb} psi  
            AOF = {self.aof} bbl/d 
                """
            
    def pwf_to_flow(self, pwf):
        pwf = np.atleast_1d(pwf)
        q = self._pwf_to_flow(pwf)
        return q
    
    def flow_to_pwf(self, q):
        q = np.atleast_1d(q)
        pwf = self._flow_to_pwf(q)
        return pwf
    
    def flow_to_dd(self, q):
        q = np.atleast_1d(q)
        pwf = self._flow_to_pwf(q)
        dd = self.pr - pwf
        return dd
    
    def dd_to_flow(self, dd):
        dd = np.atleast_1d(dd)
        pwf = self.pr - dd
        q = self._pwf_to_flow(pwf)
        return q
        
    def create_curve(self,n=10):
        c,_ = oil_inflow_curve(self.pr,self.j,self.pb,n=n)
        return c
    
    def plot(self,ax=None,n=10,**kwargs):
        
        pb = kwargs.pop('pb',True)
        q = kwargs.pop('flow',None)
        pwf = kwargs.pop('pwf',None)
        dd = kwargs.pop('dd',None)
        
        def_kw = {
        'color': 'darkgreen',
        'linestyle':'-',
        'linewidth': 2
        }    
        for (k,v) in def_kw.items():
            if k not in kwargs:
                kwargs[k]=v
                
        oax=ax or plt.gca()
        curve = self.create_curve(n)
        oax.plot(curve['q'],curve['p'],**kwargs)
        oax.set_xlabel("Flow Rate [bbl/d]")
        oax.set_ylabel("Pwf [psi]")
        
        if pb==True:
            qpb = self.pwf_to_flow(self.pb)
            oax.scatter(qpb,self.pb, color='red', label='Bubble Point')
        
        if q is not None:
            q = np.atleast_1d(q)
            p = self._flow_to_pwf(q)
            oax.scatter(q,p, color='springgreen', s=70)
            
        if pwf is not None:
            pwf = np.atleast_1d(pwf)
            flow = self._pwf_to_flow(pwf)
            oax.scatter(flow,pwf, color='lime',s=70)
            
        if dd is not None:
            dd = np.atleast_1d(dd)
            pwf_dd = self.pr - dd
            flow_dd = self._pwf_to_flow(pwf_dd)
            oax.scatter(flow_dd,pwf_dd, color='lime',s=70)
        
