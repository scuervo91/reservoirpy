import numpy as np

#####################################################################################
############################ OIL CORRELATIONS #######################################

# Bubble point Correlations

#Correction by Non-Hydrocarbon gases

#Correction by N2
def n2_correction(api=None,temp=None,y=0):
    if y==0:
        cn2 = 1
    else:
        cn2 = 1.0 + ((-2.65e-4 * api + 5.5e-3)*temp + (0.0931*api - 0.8295))*y + ((1.954e-11 * np.power(api,4.699)*temp)+(0.027*api - 2.366))*np.power(y,2)
    return cn2 

def co2_correction(y=0,temp=None):
    if y==0:
        cco2 = 1
    else:
        cco2 = 1.0 - 693.8*y*np.power(temp,-1.553)
    return cco2 

def h2s_correction(api=None,y=0,temp=None):
    if y==0:
        ch2s = 1
    else:
        ch2s = 1.0 - (0.9035 + 0.0015*api)*y + 0.019*(45-api)*np.power(y,2)
    return ch2s 

#Standing
def pb_standing(rs=None,temp=None,sg_gas=None,api=None,**kwargs):
    """
    Estimate the bubble point pressure using Standing Correlation

    Input: 
        rs -> Solution Gas Ratio [scf/bbl]
        temp -> Temperature [F]
        sg_gas -> Specific Gravity gas (air=1)
        api -> Oil API gravity [API]
    
    Return:
        pb -> Bubble Point Pressure

    Source: Correlaciones Numericas PVT - Carlos Banzer
    """
    rs = np.atleast_1d(rs)
    assert isinstance(rs,(np.ndarray))
    
    temp = np.atleast_1d(temp)
    assert isinstance(temp,(np.ndarray))
    
    sg_gas = np.atleast_1d(sg_gas)
    assert isinstance(sg_gas,(np.ndarray))
    
    api = np.atleast_1d(api)
    assert isinstance(api,(np.ndarray))
     
    #Corrections for non Hydrocarbon gases
    y_n2 = kwargs.pop('y_n2',0)
    y_co2 = kwargs.pop('y_co2',0)
    y_h2s = kwargs.pop('y_h2s',0)

    cn2 = n2_correction(y=y_n2,api=api,temp=temp)
    cco2 = co2_correction(y=y_co2,temp=temp)
    ch2s = h2s_correction(y=y_h2s,api=api,temp=temp)

    
    f = np.power(rs/sg_gas,0.83) * np.power(10,0.00091*temp - 0.0125*api)
    pb = 18.2 * (f - 1.4)
    pb_corrected = pb * cn2 * cco2 * ch2s
    return pb_corrected
