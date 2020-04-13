import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from datetime import date

############################################################################################
# Forecast Function
def forecast_curve(RangeTime,qi,di,ti,b):
  """
  Estimate a Forecast curve given Decline curve parameters

  Attributes:
    RangeTime: Range of dates to estimate the Forecast Curve-> Timestamp Series
    qi:        Initial flow rate -> Number
    di:        Decline rate in fraction and positive-> Number
    ti:        Date of the initial flow Rate-> Timestamp
    b:         Arp's Coefficient. 0<=b<=1  -> Number 

    Return -> Three-Column DataFrame: 
              -Column 'time' Timestamp Series 
              -Column 'curve' Forecast values Series
              -Column 'cum' cummulative flow rate

  """
  ##Convert dates to number for apply regression methods
  days_number = RangeTime.apply(lambda x: x.toordinal()) 
  ti_day = ti.toordinal()

  #Estimate the difference in days between the dates to forecast and Initial Ti                                
  day_diff = days_number-ti_day                          

  if b == 0:
    q = qi*np.exp(-(di/365)*day_diff) 
  elif (b>0)&(b<=1):
    q = qi/np.power(1+b*(di/365)*day_diff,1/b)

  diff_period = np.append(np.array([1]),np.diff(days_number))
  diff_q = diff_period * q 
  cum = diff_q.cumsum()
  forecast = pd.DataFrame({'time':RangeTime,'rate':q, 'cum':cum})
  forecast = forecast.set_index('time')
  Np = forecast.iloc[-1,-1]
  
  return forecast, Np

def forecast_econlimit(t,qt,qi,di,ti,b, fr):
  """
  Estimate a Forecast curve until a given Economic limit rate and Decline curve parameters

  Attributes:
    t:         Initial date to start forecast
    qt: Economic limit rate -> Number
    qi:        Initial flow rate -> Number
    di:        Decline rate in fraction and positive-> Number
    ti:        Date of the initial flow Rate-> Timestamp
    b:         Arp's Coefficient. 0<=b<=1  -> Number 

  Return -> Three-Column DataFrame: 
            -Column 'time' Timestamp Series 
            -Column 'curve' Forecast values Series
            -Column 'cum' cummulative flow rate
  """
  # Estimate the time at economic limit
  if b == 0:
    date_until = pd.Timestamp.fromordinal(int(np.log(qi / qt) * (1/(di/365))) + ti.toordinal())
  elif (b > 0) & (b <= 1):
    date_until = pd.Timestamp.fromordinal(int((np.power(qi / qt, b) - 1)/(b * (di/365))) + ti.toordinal())
  else:
    raise ValueError('b must be between 0 and 1')

  TimeRange = pd.Series(pd.date_range(start=t, end=date_until, freq=fr))

  f, Np = forecast_curve(TimeRange,qi,di,ti,b)

  return f, Np

######################################################################
#Create Declination Object 

class declination:
  """
  Decline curve object for Oil and Gas Forecasting
  
  Attributes:
    Qi: Initial Flow Rate in bbl/d: Number
    Di: Decline rate. Number must be positive and written in fraction: Number
    Ti: Date if Initial flow Rate (Qi): Timestamp
    b:  Arps Coefficient: 0<=b<=1  
  
   """
  def __init__(self, **kwargs):
    self.qi = kwargs.pop('qi',None)
    self.di = kwargs.pop('di',None)
    self.b = kwargs.pop('b',0)
    self.ti = kwargs.pop('ti',None)

#####################################################
############## Properties ###########################

  @property
  def qi(self):
    return self._qi

  @qi.setter
  def qi(self,value):
    assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be number'
    self._qi = value

  @property
  def di(self):
    return self._di

  @di.setter
  def di(self,value):
    assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be number'
    self._di = value

  @property
  def b(self):
    return self._b

  @b.setter
  def b(self,value):
    assert isinstance(value,(int,float,np.ndarray,type(None))), f'{type(value)} not accepted. Name must be number'
    assert value >= 0 and value <= 1
    self._b = value

  @property
  def ti(self):
    return self._ti

  @ti.setter
  def ti(self,value):
    assert isinstance(value,(date,type(None))), f'{type(value)} not accepted. Name must be date'
    self._ti = value

  @property
  def kind(self):
    if self._b == 0:
      self._kind='Exponential'
    elif self._b == 1:
      self._kind = 'Harmonic'
    elif (self._b<1)&(self._b>0):
      self._kind = 'Hyperbolic'
    return self._kind
    
  def __str__(self):
    return '{self.kind} Declination \n Ti: {self.ti} \n Qi: {self.qi} bbl/d \n Rate: {self.di} Annually \n b: {self.b}'.format(self=self)
  
  def __repr__(self):
    return '{self.kind} Declination \n Ti: {self.ti} \n Qi: {self.qi} bbl/d \n Rate: {self.di} Annually \n b: {self.b}'.format(self=self)

  def forecast(self,start_date, end_date, fq='M'):
    """
    Forecast curve from the declination object.

    Parameters: 
      start_date: datetime
      end_date: datetime 

    Return: 
      f: DataFrame with t column and curve column
      np: Cummulative production

    """
    assert isinstance(start_date,date), 'start_date must be date'
    assert isinstance(end_date,date), 'send_date must be date'

    TimeRange = pd.Series(pd.date_range(start=start_date, end=end_date, freq=fq, closed=None))
    f, Np = forecast_curve(TimeRange,self.qi,self.di,self.ti,self.b)
    return f, Np

  def forecast_econ(self, t, qt, fq='M'):
    """
    Estimate a Forecast curve until a given Economic limit rate and Decline curve parameters

  Attributes:
    t:         Initial date to start forecast
    qt:        Economic limit rate -> Number
    fq:       Forecast frequency-> default 'M' Monthly

  Return -> Three-Column DataFrame: 
            -Column 'time' Timestamp Series 
            -Column 'curve' Forecast values Series
            -Column 'cum' cummulative flow rate

    """
    f, Np = forecast_econlimit(t,qt,self.qi,self.di,self.ti,self.b, fr=fq)

    return f, Np


################################################################################
#Decline Fit
def decline_fit(RangeTime,FlowRate,b=None, ad=True):
  """
  Estimate the declination parameters of a time series of production daily rate
  as a Decline Curve defined by Arps

     Attributes:
    RangeTime: Range of dates to estimate the declination-> Timestamp Series
    FlowRate:  Production rate in daily basis-> Series
    b:         Arp's Coefficient. 0<=b<=1  -> If None b parameter is also fitted
                                           -> if  (b>=0)&(b<=1) b is not fitted but fixed
                                           -> Default: None
    ad:        apply anomally detection    ->  Bool, Default: True
              

    Return -> q -> 1D Numpy array with the Flow rate
  """
  if ad == True:
    tnum = RangeTime.apply(lambda x: x.toordinal()) 
    lnq = np.log(FlowRate)
    slp = -np.diff(lnq) / np.diff(tnum)
    slp = np.append(slp[0],slp)
    mu = slp.mean()
    sig=slp.std()
    RangeTime = RangeTime[np.abs(slp)<mu+2*sig]
    FlowRate = FlowRate[np.abs(slp)<mu+2*sig]

  if b is None:
    
    def decline_function(RangeTime,qi,di,b):
      """
      Estimate the flow rate given the decline curve parameters assuming the Ti 
      to the first value on RangeTime. 
      
      This function is intended to be used with scipy.optimize.curve_fit function to
      create the cost function to fit the decline curve parameters to a given 
      production data. 

      Attributes:
        RangeTime: Range of dates to estimate the Forecast Curve-> Timestamp Series
        qi:        Initial flow rate -> Number
        di:        Decline rate in fraction and positive-> Number
        b:         Arp's Coefficient. 0<=b<=1  -> Number 

        Return -> q -> 1D Numpy array with the Flow rate: 
      """
      days_number = RangeTime.apply(lambda x: x.toordinal())
      ti_day = RangeTime[0].toordinal() 
      day_diff = days_number-ti_day 

      if b == 0:
        q = qi*np.exp(-(di/365)*day_diff) 
      elif (b>0) & (b<=1):
        q = qi/np.power(1+b*(di/365)*day_diff,1/b)
      
      return q
    
    popt, pcov = curve_fit(decline_function, RangeTime, FlowRate, bounds=(0, [np.inf, np.inf, 1]))
    dec = declination(qi=popt[0], di=popt[1], ti=RangeTime[0], b=popt[2])
  elif (b >= 0) & (b <= 1):
    
    def decline_function(RangeTime,qi,di):
      """
      Estimate the flow rate given the decline curve parameters assuming the Ti 
      to the first value on RangeTime. 
      
      This function is intended to be used with scipy.optimize.curve_fit function to
      create the cost function to fit the decline curve parameters to a given 
      production data. 

      Attributes:
        RangeTime: Range of dates to estimate the Forecast Curve-> Timestamp Series
        qi:        Initial flow rate -> Number
        di:        Decline rate in fraction and positive-> Number
        b:         Arp's Coefficient. 0<=b<=1  -> Number 

        Return -> declination object
      """
      days_number = RangeTime.apply(lambda x: x.toordinal())
      ti_day  = RangeTime[0].toordinal() 
      day_diff = days_number-ti_day 
      b
      if b == 0:
        q = qi*np.exp(-(di/365)*day_diff) 
      elif (b>0)&(b<=1):
        q = qi/np.power(1+b*(di/365)*day_diff,1/b)
      
      return q 
    
    popt, pcov = curve_fit(decline_function, RangeTime, FlowRate, bounds=(0, [np.inf, np.inf]))
    dec = declination(qi=popt[0], di=popt[1], ti=RangeTime[0], b=b)
  
  return dec