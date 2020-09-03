import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from ...pvtpy.black_oil import pvt, gas



def oil_j(mu=None,bo=None,h=None,k=None,kh=None,re=1490,rw=0.58,s=0):
    """
    Estimate the Oil productivity index with the pseudosteady-state flow regime equation in bbl/d*psi
    """
    if kh is None:
        kh = k*h
    j=(0.00708*kh)/(mu*bo*(np.log(re/rw)-0.75+s))
    return j

def gas_j(h=None,k=None,kh=None,re=1490,rw=0.58,s=0, temp=None):
    """
    Estimate the Gas productivity index with the pseudosteady-state flow regime equation in Kscfd/d*psi^2*cP
    """
    #If kh is not given estimate
    if kh is None:
        kh = k*h
    
    #Convert temperature from °F to °R
    temp += 460
    j = kh/(1422*temp*(np.log(re/rw)-0.75+s))
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

def gas_inflow_curve(pr,j,gas_pvt,n=10):
    """
    Estimate the Pressure and Flow curves given a Productivity Index parameters

    Attributes:
        pr:     Reservoir Pressure [psi]
        j:      Productivity index [bbl/d*psi]
        pb:     PVT
        n:      Number of steps in the return DataFrame
    
    Return:
        pi:     DataFrame with two columns Pressure [p] and Flow [q]-> Pandas DataFrame
        aof:    Absolute open flow -> Number
    """
    assert isinstance(gas_pvt, pvt), "The pvt must be of PVT class"
    assert all(i in gas_pvt.columns for i in ['z','mug'])

    #Create pressure Range
    p = np.linspace(0,pr,n)

    #Estimate pvt interpolated
    pvt_interpolated = gas_pvt.interpolate(p)

    #Calculate 2p/(mug*z) Pseudopressure
    pvt_interpolated['pspressure'] = (2*pvt_interpolated.index) / (pvt_interpolated['z']*pvt_interpolated['mug'])

    q=np.zeros(p.shape)

    #Calculate flow
    mp_pr = np.trapz(pvt_interpolated['pspressure'],pvt_interpolated.index)
    for i in range(1,pvt_interpolated.shape[0]):
        mp_pwf = np.trapz(pvt_interpolated['pspressure'].iloc[:i],pvt_interpolated.index[:i])
        q[i-1] = j * (mp_pr - mp_pwf)

    data = pd.DataFrame({'p':p,'q':q})
    aof = data['q'].max()

    return data, aof




class oil_inflow:

    def __init__(self,**kwargs):
        self.pr = kwargs.pop('pr',None)
        self.j = kwargs.pop('j',None)
        self.pb = kwargs.pop('pb',None) 
        self.n = kwargs.pop('n',20) 


#####################################################
############## Properties ###########################

    @property
    def pr(self):
        return self._pr

    @pr.setter
    def pr(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be number'
        self._pr = value

    @property
    def j(self):
        return self._j

    @j.setter
    def j(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be number'
        self._j = value

    @property
    def pb(self):
        return self._pb

    @pb.setter
    def pb(self,value):
        if value is not None:
            assert isinstance(value,(int,float,np.ndarray)), f'{type(value)} not accepted. Name must be number'
        else:
            value = self.pr
        self._pb = value

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self,value):
        assert isinstance(value,int), f'{type(value)} not accepted. Name must be int'
        self._n = value

    @property
    def aof(self):
        _,_aof = oil_inflow_curve(self._pr,self._j,self.pb,n=self._n)
        return _aof

    @property
    def df(self):
        _df,_ = oil_inflow_curve(self._pr,self._j,self.pb,n=self._n)
        return _df

#####################################################
############## repr ###########################

    def __str__(self):
        return f"""Oil Inflow: 
            Reservoir Pressure: {self.pr} psi 
            Productivity Index: {round(self.j,2)} bbl/d*psi 
            Bubble Point Pressure: {self.pb} psi  
            AOF = {round(self.aof,2)} bbl/d 
                """
  
    def __repr__(self):
        return f"""Oil Inflow: 
            Reservoir Pressure: {self.pr} psi 
            Productivity Index: {round(self.j,2)} bbl/d*psi 
            Bubble Point Pressure: {self.pb} psi  
            AOF = {round(self.aof,2)} bbl/d 
                """


#####################################################
############## methods ###########################
            
    def pwf_to_flow(self, pwf):
        _df = self.df
        _pwf_to_flow = interp1d(_df['p'],_df['q'])
        pwf = np.atleast_1d(pwf)
        q = _pwf_to_flow(pwf)
        return q
    
    def flow_to_pwf(self, q):
        _df = self.df
        _flow_to_pwf = interp1d(_df['q'],_df['p'])
        q = np.atleast_1d(q)
        pwf = _flow_to_pwf(q)
        return pwf
    
    def flow_to_dd(self, q):
        _df = self.df
        _flow_to_pwf = interp1d(_df['q'],_df['p'])        
        q = np.atleast_1d(q)
        pwf = _flow_to_pwf(q)
        dd = self._pr - pwf
        return dd
    
    def dd_to_flow(self, dd):
        _df = self.df
        _pwf_to_flow = interp1d(_df['p'],_df['q'])
        dd = np.atleast_1d(dd)
        pwf = self._pr - dd
        q = _pwf_to_flow(pwf)
        return q

    
    def plot(self,ax=None,**kwargs):
        
        pb = kwargs.pop('pb',True)
        flow = kwargs.pop('flow',None)
        pwf = kwargs.pop('pwf',None)
        dd = kwargs.pop('dd',None)
        legend = kwargs.pop('legend',True)
        
        def_kw = {
        'color': 'darkgreen',
        'linestyle':'-',
        'linewidth': 2
        }    
        for (k,v) in def_kw.items():
            if k not in kwargs:
                kwargs[k]=v
                
        oax = ax or plt.gca()
        _df = self.df
        _flow_to_pwf = interp1d(_df['q'],_df['p']) 
        _pwf_to_flow = interp1d(_df['p'],_df['q'])
        oax.plot(_df['q'],_df['p'],**kwargs)
        oax.set_xlabel("Flow Rate [bbl/d]")
        oax.set_ylabel("Pwf [psi]")
        
        if pb==True:
            qpb = _pwf_to_flow(self._pb)
            oax.scatter(qpb,self._pb, color='red', label='Bubble Point')
        
        if flow is not None:
            flow = np.atleast_1d(flow)
            pwf_f = _flow_to_pwf(flow)
            oax.scatter(flow,pwf_f, color='springgreen', s=70,label=f" Flow {flow} bbl/d \n Pwf{np.round(pwf_f,decimals=1)} psi")
            
        if pwf is not None:
            pwf = np.atleast_1d(pwf)
            flow_p = _pwf_to_flow(pwf)
            oax.scatter(flow_p,pwf, color='lime',s=70, label=f" Pwf {pwf} psi \n Flow{np.round(flow_p,decimals=1)} bbl/d")
            
        if dd is not None:
            dd = np.atleast_1d(dd)
            pwf_dd = self._pr - dd
            flow_dd = _pwf_to_flow(pwf_dd)
            oax.scatter(flow_dd,pwf_dd, color='lime',s=70,label=f" DD {dd} psi \n Pwf{np.round(pwf_dd,decimals=1)} psi \n Flow{np.round(flow_dd,decimals=1)} bbl/d")

        if legend:
            oax.legend()

class gas_inflow:

    def __init__(self,**kwargs):
        self.pr = kwargs.pop('pr',None)
        self.j = kwargs.pop('j',None)
        self.gas = kwargs.pop('gas',None) 
        self.n = kwargs.pop('n',20) 


#####################################################
############## Properties ###########################

    @property
    def pr(self):
        return self._pr

    @pr.setter
    def pr(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be number'
        self._pr = value

    @property
    def j(self):
        return self._j

    @j.setter
    def j(self,value):
        assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be number'
        self._j = value

    @property
    def gas(self):
        return self._gas

    @gas.setter
    def gas(self,value):
        assert isinstance(value,gas), f'{type(value)} not accepted. Name must be gas type'
        self._gas = value

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self,value):
        assert isinstance(value,int), f'{type(value)} not accepted. Name must be int'
        self._n = value

    @property
    def aof(self):
        _,_aof = gas_inflow_curve(self.pr,self.j,self.gas.pvt,n=self._n)
        return _aof

    @property
    def df(self):
        _df,_ = gas_inflow_curve(self.pr,self.j,self.gas.pvt,n=self._n)
        return _df



#####################################################
############## repr ###########################

    def __str__(self):
        return f"""Gas Inflow: 
            Reservoir Pressure: {self.pr} psi 
            Productivity Index: {round(self.j,2)} Kscf/d*psi^2*cP
            AOF = {round(self.aof,2)} Kscf/d 
                """
  
    def __repr__(self):
        return f"""Gas Inflow: 
            Reservoir Pressure: {self.pr} psi 
            Productivity Index: {round(self.j,2)} Kscf/d*psi^2*cP
            AOF = {round(self.aof,2)} Kscf/d 
                """


#####################################################
############## methods ###########################
            
    def pwf_to_flow(self, pwf):
        _df = self.df
        _pwf_to_flow = interp1d(_df['p'],_df['q'])
        pwf = np.atleast_1d(pwf)
        q = _pwf_to_flow(pwf)
        return q
    
    def flow_to_pwf(self, q):
        _df = self.df
        _flow_to_pwf = interp1d(_df['q'],_df['p'])
        q = np.atleast_1d(q)
        pwf = _flow_to_pwf(q)
        return pwf
    
    def flow_to_dd(self, q):
        _df = self.df
        _flow_to_pwf = interp1d(_df['q'],_df['p'])        
        q = np.atleast_1d(q)
        pwf = _flow_to_pwf(q)
        dd = self._pr - pwf
        return dd
    
    def dd_to_flow(self, dd):
        _df = self.df
        _pwf_to_flow = interp1d(_df['p'],_df['q'])
        dd = np.atleast_1d(dd)
        pwf = self._pr - dd
        q = _pwf_to_flow(pwf)
        return q

    
    def plot(self,ax=None,**kwargs):
        
        flow = kwargs.pop('flow',None)
        pwf = kwargs.pop('pwf',None)
        dd = kwargs.pop('dd',None)
        legend = kwargs.pop('legend',True)
        
        def_kw = {
        'color': 'darkred',
        'linestyle':'-',
        'linewidth': 2
        }    
        for (k,v) in def_kw.items():
            if k not in kwargs:
                kwargs[k]=v
                
        oax = ax or plt.gca()
        _df = self.df
        _flow_to_pwf = interp1d(_df['q'],_df['p']) 
        _pwf_to_flow = interp1d(_df['p'],_df['q'])
        oax.plot(_df['q'],_df['p'],**kwargs)
        oax.set_xlabel("Flow Rate [kscf/d]")
        oax.set_ylabel("Pwf [psi]")
        
      
        if flow is not None:
            flow = np.atleast_1d(flow)
            pwf_f = _flow_to_pwf(flow)
            oax.scatter(flow,pwf_f, color='springgreen', s=70,label=f" Flow {flow} kscfd \n Pwf{np.round(pwf_f,decimals=1)} psi")
            
        if pwf is not None:
            pwf = np.atleast_1d(pwf)
            flow_p = _pwf_to_flow(pwf)
            oax.scatter(flow_p,pwf, color='lime',s=70, label=f" Pwf {pwf} psi \n Flow{np.round(flow_p,decimals=1)} kscfd")
            
        if dd is not None:
            dd = np.atleast_1d(dd)
            pwf_dd = self._pr - dd
            flow_dd = _pwf_to_flow(pwf_dd)
            oax.scatter(flow_dd,pwf_dd, color='lime',s=70,label=f" DD {dd} psi \n Pwf{np.round(pwf_dd,decimals=1)} psi \n Flow{np.round(flow_dd,decimals=1)} bbl/d")

        if legend:
            oax.legend()