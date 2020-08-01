import numpy as np 
import pandas as pd 
from scipy.interpolate import interp1d

def vshale_gr(gr_curve,gr_sand,gr_shale,type='linear'):
    gr_curve=np.atleast_1d(gr_curve)
    gr_sand=np.atleast_1d(gr_sand)
    gr_shale=np.atleast_1d(gr_shale)
    igr=(gr_curve-gr_sand)/(gr_shale-gr_sand)
    igr[igr < 0.0] = 0.0
    igr[igr > 1.0] = 1.0
    #https://www.geoloil.com/VshModels.php
    if type == 'linear':
        vsh = igr
    elif type == 'clavier':
        vsh = 1.7 - np.sqrt(3.38 - np.power(igr+0.7,2))
    elif type == 'stieber':
        vsh = igr/(3-2*igr)
    elif type == 'larionov_old':
        vsh = 0.33 * (np.power(2,2*igr)-1)
    elif type == 'larionov_tertiary':
        vsh = 0.083 * (np.power(2,3.7*igr)-1)
    else:
        raise ValueError(f'method especified [ {type} ] does not exist')

    return vsh

def vshale_dn(rho_curve, ntr_curve, rho_ma=2.65, rho_f=1.0, hi_shl=0.46,rho_shl=2.43):
    rho_curve= np.atleast_1d(rho_curve)
    ntr_curve= np.atleast_1d(ntr_curve)
    rho_ma = np.atleast_1d(rho_ma)
    rho_f = np.atleast_1d(rho_f)
    hi_shl = np.atleast_1d(hi_shl)
    rho_shl = np.atleast_1d(rho_shl)
    vsh = (rho_curve - rho_ma + ntr_curve*(rho_ma-rho_f))/(rho_shl - rho_ma + hi_shl*(rho_ma-rho_f))
    vsh[vsh < 0.0] = 0.0
    vsh[vsh > 1.0] = 1.0
    return vsh
    

def phi_rho(rho_curve,rho_ma=2.65,rho_f=1.0):
    rho_curve=np.atleast_1d(rho_curve)
    rho_ma=np.atleast_1d(rho_ma)
    rho_f=np.atleast_1d(rho_f)
    phi_rho_curve=(rho_ma-rho_curve)/(rho_ma-rho_f)
    phi_rho_curve[phi_rho_curve < 0.0] = 0.0
    phi_rho_curve[phi_rho_curve > 1.0] = 1.0
    return phi_rho_curve

def phie(phi_curve,vsh_curve):
    phi_curve=np.atleast_1d(phi_curve)
    vsh_curve=np.atleast_1d(vsh_curve)
    phie_curve=phi_curve*(1 -vsh_curve)
    phie_curve[phie_curve < 0.0] = 0.0
    phie_curve[phie_curve > 0.3] = 0.3
    return phie_curve

def phia(phi_rho_curve, ntr_curve, method='geometric'):
    phi_rho_curve = np.atleast_1d(phi_rho_curve)
    ntr_curve = np.atleast_1d(ntr_curve)
    c = np.transpose(np.vstack((phi_rho_curve,ntr_curve)))
    if method == 'mean':
        phia_curve = np.mean(c,axis=1)
    elif method== 'geometric':
        phia_curve = np.power(((np.power(phi_rho_curve,2)+np.power(ntr_curve,2))/2),0.5)
    return phia_curve
        
def facies_dnp(rho_curve, ntr_curve,pef_curve,**kw):
    rho_curve = np.atleast_1d(rho_curve)
    ntr_curve = np.atleast_1d(ntr_curve)
    pef_curve = np.atleast_1d(pef_curve)
    phi_rho_curve = phi_rho(rho_curve,**kw)
    phia_curve = phia(phi_rho_curve,ntr_curve)
    
    u = pef_curve*((rho_curve + 0.1833)/1.07)
    uma = (u - 0.398 * phia_curve)/(1-phia_curve)
    dga = (rho_curve - phia_curve)/(1-phia_curve)
    return uma, dga
    
def sw(rt_curve,phi_curve,rw,vsh_curve=None,a=0.62,m=2.15,n=2,rsh=4.0,alpha=0.3,method="archie"):
    a=np.atleast_1d(a)
    m=np.atleast_1d(m)
    n=np.atleast_1d(n)
    vsh = np.atleast_1d(vsh_curve) if vsh_curve is not None else None
    rsh=np.atleast_1d(rsh)
    alpha=np.atleast_1d(alpha)
    rt=np.atleast_1d(rt_curve)
    phi = np.atleast_1d(phi_curve)
    rw=np.atleast_1d(rw)
    if method == "archie":
        sw_curve=np.power(((a*rw)/(rt*np.power(phi,m))),1/n)
    elif method == "smdx": #https://www.spec2000.net/14-sws.htm
        C=((1-vsh)*a*rw)/np.power(phi,m)
        D=C*vsh/(2 * rsh)
        E=C/rt
        sw_curve=np.power(np.sqrt(D**2 + E) - D, 2/n)
    elif method == "indo":
        #https://geoloil.com/Indonesia_SW.php
        #A=np.sqrt(1 /rt)
        #B=(np.power(vsh,(1 -(vsh/2)))/np.sqrt(rsh))
        #C=np.sqrt(np.power(phi,m)/(a*rw))
        #sw_curve=np.power((A/(B+C)),2/n)

        #http://nafta.wiki/display/GLOSSARY/Indonesia+Model+%28Poupon-Leveaux%29+@model
        A_inv = 1 + np.sqrt((np.power(vsh,2-vsh)*rw)/(phi*rsh))
        A = 1/A_inv
        sw_curve = np.power((A*rw)/(np.power(phi,m)*rt),1/n)
    elif method == "fertl":
        A=np.power(phi,-m/2)
        B=(a*rw)/rt
        C=((alpha*vsh)/2)**2
        sw_curve=A*((np.sqrt(B+C))-np.sqrt(C))
    sw_curve[sw_curve < 0.0] = 0.0
    sw_curve[sw_curve > 1.0] = 1.0
    
    return sw_curve

def depth_temperature(depth, surface_temperature=77 ,gradient=1):
    depth = np.atleast_1d(depth)
    
    t = (gradient/100) * depth + surface_temperature 
    return t

def rw_temp_convert(rw,t1,t2, temp_unit='f'):
    rw = np.atleast_1d(rw)
    t1 = np.atleast_1d(t1)
    t2 = np.atleast_1d(t2)
    
    if temp_unit=='f':
        c = np.array([6.77])
    else:
        c = np.array([21.5])
    
    rw2 = rw*((t1 + c)/(t2 + c))
    return rw2

def rw(temp, salinity,temp_unit='f'):
    """
    Tc = 60.0       # Temperature (F)
    Wse =  60000.0  # Salinity (ppm)
    """
    
    # 1) Convert from Celcius to Farenheit
    if temp_unit=='c':
        tf = 1.8*temp + 32.0
    else:
        tf=temp
    
    # 2) Calculate Resistivity in Ohm meters
    rw = np.power((400000.0/(tf*salinity)),0.88)
    
    return rw



def rw2(T, Wse, verbose=False,Celcius=False):    
    """
    Uses method from textbook "Petrophysics" by Djebbar and Donaldson.
    
    Supposedly more accurate than the approach used in calc_Rw because Hx and 
    RwT account for non-linearlity in resistivity as a function of salinity.
    
    
    #  Input Parameters
    # -------------------
    
    Tc = 60.0       # Temperature (Celcius)
    Wse = 60000.0   # Salinity (ppm)
    """
    
    #  Calculations:
    
    # 1) Convert from Celcius to Farenheit
    if Celcius==True:
        Tf = 9.0/5.0*T + 32.0
    else:
        Tf=T
    
    # 2) Calculate reference water resistivity @ 75 Degrees Farenheit
    Rw75 = 1.0/(2.74*10**-4 * Wse**0.955) + 0.0123
    
    # 3) Calculate nonlinear correction factors
    Xh = 10.0**(-0.3404*np.log10(Rw75) + 0.6414)
    
    # 4) Calculate Water Resistivity at Temperature T1.  Output is Ohm-m
    Rw = Rw75 * (75.0 + Xh)/(Tf + Xh)
    
    if verbose == True:
        print(" ")
        print("T (f):     %10.2f" % (Tf))
        print("Wse (ppm): %10.2f" % (Wse))
        print("Rw (Ohm*m):%10.5f" % (Rw))
    
    return Rw

def perm(phie_curve,swir,fluid='oil',author='timur'):
    phie=np.atleast_1d(phie_curve)
    swir=np.atleast_1d(swir)
    if author=='timur':
        Dperm=4.4
        Eperm=2.0

        if fluid=="oil":
            Cperm=3400
        elif fluid=="gas":
            Cperm=340
        k=Cperm * np.power(phie,Dperm) / np.power(swir,Eperm)
        
    if author=='morris':
        Dperm=6.0
        Eperm=2.0

        if fluid=="oil":
            Cperm=65000
        elif fluid=="gas":
            Cperm=6500
        k=Cperm * np.power(phie,Dperm) / np.power(swir,Eperm) 
    
    if author=='coates':
        k=(10.0*phie)**4 *((1.0-swir)/swir)
        
    return k


def flow_capacity(height,perm_curve,pay_curve):
    perm_curve=np.nan_to_num(np.atleast_1d(perm_curve))
    pay_curve=np.nan_to_num(np.atleast_1d(pay_curve))
    height=np.atleast_1d(height)
    
    kh=height*perm_curve*pay_curve
    khcum=np.cumsum(kh)
    kht=np.sum(kh)
    khnorm=1-(khcum/kht)
    return kh, khnorm

def sw_pnn(phie,vsh, sigma,sighy,sigsh,ws=None,sigw=None,sigmam=None):
    """
    https://www.spec2000.net/14-swtdt.htm
    PHIe = effective porosity (fractional)
    SIGMA = TDT capture cross section log reading (capture units)
    SIGMAM = capture cross section matrix value (capture units)
    SIGW = capture cross section for water (capture units)
    SIGHY = capture cross section for hydrocarbons (capture units)
    SIGSH = capture cross section for shale (capture units)
    SWtdt = water saturation from TDT (fractional)
    Vsh = shale volume (fractional)
    WS = water salinity (ppm NaCl)
    """
    phie=np.atleast_1d(phie)
    vsh=np.atleast_1d(vsh)
    sigma=np.atleast_1d(sigma)
    sighy=np.atleast_1d(sighy)
    sigsh=np.atleast_1d(sigsh)

    if sigw is None:
        sigw = 22 + 0.000404*ws 
    if sigmam is None:
        sigmam = (sigma - phie*sigw)/(1-phie)

    _a = sigma - sigmam 
    _b = phie*(sighy - sigmam)
    _c = vsh*(sigsh - sigmam)
    _d = phie*(sigw - sighy)
    sw = (_a - _b - _c)/(_d)
    sw[phie<=0] = 1.0

    return sw

def rw_from_sp(rmf=None, rmf_temp=75, res_temp=None, ssp=None,temp_unit='f', rw_75=False):
    """
    Estimate water resistivity from SP log
    """

    #https://www.spec2000.net/05-7rwsp.htm
    # 1) Convert from Celcius to Farenheit
    if temp_unit=='c':
        rmf_temp = 1.8*rmf_temp + 32.0
        res_temp = 1.8*res_temp + 32.0

    # Convert rmf @ rmf_temp to res_temp
    rmf_res = rw_temp_convert(rmf,rmf_temp,res_temp, temp_unit='f')

    #Estimate Mud filtrate equivalent resistivity
    if rmf_res > 0.1:
        rmfe = rmf_res*0.85
    else:
        rmfe = (146*rmf_res-5)/(337*rmf_res + 77)
    #Ksp
    ksp = 60 +0.122*res_temp

    #RSP
    rsp = np.power(10,-ssp/ksp)
    # Rwe. Equivalent Water Resistivity
    rwe = rmfe / rsp

    # Estimate Water resistivity from Equivalent water resistivity
    if rwe > 0.12:
        rw = np.power(10,0.69*rwe-0.24) - 0.58
    else:
        rw = (77*rwe + 5)/(146 - 337*rwe)

    if rw_75:
        rw = rw_temp_convert(rw,res_temp,75, temp_unit='f')

    return rw



    


