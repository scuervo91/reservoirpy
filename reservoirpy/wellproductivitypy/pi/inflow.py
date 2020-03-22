import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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
        self.flow = interp1d(self.curve['p'],self.curve['q'])
        self.pwf = interp1d(self.curve['q'],self.curve['p'])
    
    def create_curve(self,n=10):
        c,_ = oil_inflow_curve(self.pr,self.j,self.pb,n=n)
        return c
