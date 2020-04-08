import numpy as np
import pandas as pd
from numpy.polynomial.polynomial import polyval
from scipy.interpolate import interp1d

#####################################################################################
#####################################################################################
############################ OIL CORRELATIONS #######################################

def api_to_sg(api):
    api = np.atleast_1d(api)
    sg = 141.5/(131.5+api)
    return sg

def sg_to_api(sg):
    sg = np.atleast_1d(sg)
    api = (141.5/sg)-131.5
    return api

#####################################################################################
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
def pb(rs=None,temp=None,sg_gas=None,api=None,methods=None,**kwargs):
    """
    Estimate the bubble point pressure using Correlations

    Input: 
        rs -> Solution Gas Ratio [scf/bbl]
        temp -> Temperature [F]
        sg_gas -> Specific Gravity gas (air=1)
        api -> Oil API gravity [API]
        method -> List of correlation methods
                  ['standing',laster,'vazquez_beggs','glaso']
    
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

    assert isinstance(methods,list)

    #Corrections for non Hydrocarbon gases
    y_n2 = kwargs.pop('y_n2',0)
    y_co2 = kwargs.pop('y_co2',0)
    y_h2s = kwargs.pop('y_h2s',0)

    cn2 = n2_correction(y=y_n2,api=api,temp=temp)
    cco2 = co2_correction(y=y_co2,temp=temp)
    ch2s = h2s_correction(y=y_h2s,api=api,temp=temp)

    pb_dict = {}
     
    if 'standing' in methods:
        f = np.power(rs/sg_gas,0.83) * np.power(10,0.00091*temp - 0.0125*api)
        pb_standing = 18.2 * (f - 1.4)
        pb_standing_corrected = pb_standing * cn2 * cco2 * ch2s
        pb_dict['standing'] = pb_standing_corrected

    if 'laster' in methods:
        #estimate oil effective molecular weight
        mo = np.zeros(api.shape)
        mo[api<=40] = 630-10*api[api<=40]
        mo[api>40] = 73110*np.power(api[api>40],-1.562)


        #estimate system gas molar fraction
        sg_oil = api_to_sg(api) 
        yg = ((rs/379.3))/((rs/379.3)+(350*sg_oil/mo))


        pb_factor = np.zeros(yg.shape)
        pb_factor[yg<=0.6] = 0.679*np.exp(2.786*yg[yg<=0.6])-0.323
        pb_factor[yg>0.6] = 8.26 * np.power(yg[yg>0.6],3.56) +1.95 

        temp_r = temp + 459.67
        pb_laster = pb_factor * temp_r/sg_gas

        pb_laster_corrected = pb_laster * cn2 * cco2 * ch2s
        pb_dict['laster'] = pb_laster_corrected
    
    if 'vazquez_beggs' in methods:
        #Constants c1, c2, c3
        c1 = np.zeros(api.shape)
        c2 = np.zeros(api.shape)
        c3 = np.zeros(api.shape)

        c1[api<=30] = 0.0362
        c1[api>30] = 0.0178

        c2[api<=30] = 1.0937
        c2[api>30] = 1.1870

        c3[api<=30] = 25.724
        c3[api>30] = 23.931

        pb_vasquez = np.power(rs/(c1*sg_gas*np.exp((c3*api)/(temp+460))),1/c2)
        pb_vazquez_corrected = pb_vasquez * cn2 * cco2 * ch2s
        pb_dict['vazquez_beggs'] = pb_laster_corrected

    if 'glaso' in methods:
        f = np.power(rs/sg_gas,0.816)*((np.power(temp,0.172))/(np.power(api,0.989)))

        pb_glaso = np.power(10,polyval(np.log10(f),[1.7669,1.7447,-0.30218]))
        pb_glaso_corrected = pb_glaso* cn2 * cco2 * ch2s
        pb_dict['glaso'] = pb_glaso_corrected

    return pb_dict


#####################################################################################
# gas-Oil Ratio Correlations

def rs(p=None,rs=None,temp=None,sg_gas=None,api=None,pb=None,methods=None,**kwargs):
    """
    Estimate the Gas-Oil Ratio using Standing Correlation

    Input: 
        p -> Interest Pressure [psi]
        pb -> Bubble Point [psi]
        temp -> Temperature [F]
        api -> Oil API gravity [API]
        sg_gas -> Gas specifi gravity
        method -> List of correlation methods
                  ['standing',laster,'vazquez_beggs','glaso','valarde']
    
    Return:
        rs -> Gas Oil Ratio

    Source: Correlaciones Numericas PVT - Carlos Banzer
    """
    p = np.atleast_1d(p)
    assert isinstance(p,(np.ndarray))
    
    temp = np.atleast_1d(temp)
    assert isinstance(temp,(np.ndarray))
      
    api = np.atleast_1d(api)
    assert isinstance(api,(np.ndarray))

    sg_gas = np.atleast_1d(sg_gas)
    assert isinstance(sg_gas,(np.ndarray))

    pb = np.atleast_1d(pb)
    assert isinstance(pb,(np.ndarray))

    assert isinstance(methods,list)

    rs_dict = {}
    p_sat = np.zeros(p.shape)
    p_sat[p>=pb] = pb
    p_sat[p<pb] = p[p<pb]
     
    if 'standing' in methods:

        rs_standing = sg_gas * np.power(((p_sat/18.2)+1.4)*np.power(10,0.0125*api-0.00091*temp),1.2048)
        rs_dict['standing'] = rs_standing

    if 'laster' in methods:
        array_shape = p_sat * sg_gas * temp * p
        mo = np.zeros(api.shape)
        mo[api<=40] = 630-10*api[api<=40]
        mo[api>40] = 73110*np.power(api[api>40],-1.562)

        pb_factor = (p_sat*sg_gas)/(temp + 459.67)

        #estimate yg
        
        yg = np.zeros(array_shape.shape)

        yg[pb_factor<3.29] = 0.359*np.log(1.473*pb_factor[pb_factor<3.29]+0.476)
        yg[pb_factor>=3.29] = np.power(0.121*pb_factor[pb_factor>=3.29] - 0.236,0.281)

        sg_oil = api_to_sg(api)
        rs_laster = (132755*sg_oil*yg)/(mo*(1-yg))
        rs_dict['laster'] = rs_laster

    if 'vazquez_beggs' in methods:
        #Constants c1, c2, c3
        c1 = np.zeros(api.shape)
        c2 = np.zeros(api.shape)
        c3 = np.zeros(api.shape)

        c1[api<=30] = 0.0362
        c1[api>30] = 0.0178

        c2[api<=30] = 1.0937
        c2[api>30] = 1.1870

        c3[api<=30] = 25.724
        c3[api>30] = 23.931
        
        rs_vazquez = c1*sg_gas*np.power(p_sat,c2)*np.exp((c3*api)/(temp+460))
        rs_dict['vazquez_begss'] = rs_vazquez

    if 'glaso' in methods:

        f = np.power(10,2.8869-np.power(14.1811-3.3093*np.log10(p_sat),0.5))
        rs_glaso = sg_gas*np.power(f*(np.power(api,0.989)/np.power(temp,0.172)),1.2255)
        rs_dict['glaso'] = rs_glaso

    if 'valarde' in methods:
        """Method for build rs at pressures below pb by giving the rsb
        Correlation of Black Oil Properties at
        Pressures Below Bubble Point Pressure
        -A New Approach"""
        "https://wiki.pengtools.com/index.php?title=Velarde_correlation"
        A0 = 9.73e-7
        A1 = 1.672608
        A2 = 0.929870
        A3 = 0.247235
        A4 = 1.056052
        alpha_1 = A0 * np.power(sg_gas,A1)*np.power(api,A2)*np.power(temp,A3)*np.power(pb,A4)

        B0 = 0.022339
        B1 = -1.00475
        B2 = 0.337711
        B3 = 0.132795
        B4 = 0.302065
        alpha_2 = B0 * np.power(sg_gas,B1)*np.power(api,B2)*np.power(temp,B3)*np.power(pb,B4)

        C0 = 0.725167
        C1 = -1.48548
        C2 = -0.164741
        C3 = -0.09133
        C4 = 0.047094
        alpha_3 = C0 * np.power(sg_gas,C1)*np.power(api,C2)*np.power(temp,C3)*np.power(pb,C4)

        rsb = kwargs.pop('rsb',None)

        pr = p_sat/pb
        rsr = alpha_1*np.power(pr,alpha_2) + (1-alpha_1)*np.power(pr,alpha_3)
        rs_valarde = rsr*rsb
        rs_dict['valarde'] = rs_valarde

    
    rs_df = pd.DataFrame(rs_dict,index=p)

    return rs_df

#####################################################################################
#Oil Volumetric Factor

def bo(p=None,rs=None,temp=None,sg_gas=None,api=None,pb=None,methods=None,**kwargs):
    """
    Estimate the Oil Volumetric Factor using Correlations

    Input: 
        p -> Interest Pressure [psi]
        rs -> Gas Oil Ratio scf/bbl
        pb -> Bubble Point [psi]
        temp -> Temperature [F]
        api -> Oil API gravity [API]
        sg_gas -> Gas specifi gravity
        method -> List of correlation methods
                  ['standing',laster,'vazquez_beggs','glaso','valarde]
    
    Return:
        bo -> Oil Volumetric Factor

    Source: Correlaciones Numericas PVT - Carlos Banzer
    """
    p = np.atleast_1d(p)
    assert isinstance(p,(np.ndarray))

    rs = np.atleast_1d(rs)
    assert isinstance(rs,(np.ndarray))
    
    temp = np.atleast_1d(temp)
    assert isinstance(temp,(np.ndarray))
      
    api = np.atleast_1d(api)
    assert isinstance(api,(np.ndarray))

    sg_gas = np.atleast_1d(sg_gas)
    assert isinstance(sg_gas,(np.ndarray))

    pb = np.atleast_1d(pb)
    assert isinstance(pb,(np.ndarray))

    assert isinstance(methods,list)

    bo_dict = {}

    if 'standing' in methods:
        sg_oil = api_to_sg(api)
        f = rs*np.sqrt(sg_gas/sg_oil)+1.25*temp
        bo_standing = 0.9759+12e-5*np.power(f,1.2)
        bo_dict['standing'] = bo_standing

    

    if 'vazquez_beggs' in methods:
        #Constants c1, c2, c3
        c1 = np.zeros(api.shape)
        c2 = np.zeros(api.shape)
        c3 = np.zeros(api.shape)

        c1[api<=30] = 4.677e-4
        c1[api>30] = 4.670e-4

        c2[api<=30] = 1.751e-5
        c2[api>30] = 1.1e-5

        c3[api<=30] = -1.8106e-8
        c3[api>30] = 1.3370e-9

        bo_vazquez = 1+ c1*rs + c2*(temp-60)*(api/sg_gas) + c3*rs*(temp-60)*(api/sg_gas)
        bo_dict['vazquez_beggs'] = bo_vazquez

    if 'glaso' in methods:
        sg_oil = api_to_sg(api)
        f = rs*np.power(sg_gas/sg_oil,0.526) + 0.968*temp
        bo_glaso = 1 + np.power(10,-6.58511 + 2.91329* np.log10(f) - 0.27683*np.power(np.log10(f),2))
        bo_dict['glaso'] = bo_glaso

    bo_df = pd.DataFrame(bo_dict,index=p)
    return bo_df

#####################################################################################
#Oil density

def pho_oil(p=None,co=None,bo=None,rs=None,api=None,pb=None,**kwargs):
    """
    Estimate the Oil Density in lb/ft3

    Input: 
        p -> Interest Pressure [psi]
        rs -> Gas Oil Ratio scf/bbl
        pb -> Bubble Point [psi]
        co -> Isotermic oil compressibility 1/psi
        api -> Oil API gravity [API]
        Bo -> Oil Volumetric Factor

    
    Return:
        rho -> Oil Density

    Source: Correlaciones Numericas PVT - Carlos Banzer
    """
    p = np.atleast_1d(p)
    assert isinstance(p,(np.ndarray))

    rs = np.atleast_1d(rs)
    assert isinstance(rs,(np.ndarray))
    
    co = np.atleast_1d(co)
    assert isinstance(co,(np.ndarray))
      
    api = np.atleast_1d(api)
    assert isinstance(api,(np.ndarray))

    bo = np.atleast_1d(bo)
    assert isinstance(bo,(np.ndarray))

    pb = np.atleast_1d(pb)
    assert isinstance(pb,(np.ndarray))

    assert p.shape == bo.shape == rs.shape == co.shape


    rho_oil_dict = {}

    #Gas disolved specific gravity
    ygd = ((12.5+api)/50)-3.5715e-6*api*rs

    rho_oil = np.zeros(p.shape)
    p_sat = np.zeros(p.shape)
    p_sat[p>=pb] = pb
    p_sat[p<pb] = p[p<pb]

    sg_oil = api_to_sg(api)
    rho_oil[p<=pb] = (350*sg_oil+0.0764*ygd*rs[p<=pb])/(5.615*bo[p<=pb])
    
    rs_int = interp1d(p,rs)
    bo_int = interp1d(p,bo)

    rho_ob = (350*sg_oil+0.0764*ygd*rs_int(pb))/(5.615*bo_int(pb))
    rho_oil[p>pb] = rho_ob*np.exp(co[p>pb]*(pb-p))
    rho_oil_dict['density'] = rho_oil

    rho_df = pd.DataFrame(rho_oil_dict,index=p)
    return rho_df





    








