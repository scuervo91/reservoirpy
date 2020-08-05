import numpy as np
import pandas as pd 
from ...pvtpy.black_oil import pvt, gas

## Incompressible pressure drop
def potential_energy_change(
    z1=None, 
    z2=None, 
    delta_z=None,
    length=None, 
    ge=1, 
    angle=None, 
    inc=None,
    p1=0):
    """
    Δp PE accounts for the pressure change due to the weight of the column of fluid (the hydrostatic head); it
    will be zero for flow in a horizontal pipe.

    In this equation, Δz is the difference in elevation between positions 1 and 2, with z increasing upward. θ
    is defined as the angle between horizontal and the direction of flow. Thus, θ is +90° for upward, vertical
    flow, 0° for horizontal flow, and –90° for downward flow in a vertical well (Figure 7-4). For flow in a
    straight pipe of length L with flow direction θ,
    """

    # Assert height difference types
    if delta_z is None:
        if length is None:
            assert isinstance(z1,(float,int,np.ndarray)) and isinstance(z2,(float,int,np.ndarray))
            z1 = np.atleast_1d(z1)
            z2 = np.atleast_1d(z2)
            #assert z1.shape == (1,) and z2.shape == (1,)
            delta_z = z1-z2

        else:
            assert isinstance(length,(float,int,np.ndarray)) 
            length = np.atleast_1d(length)
            #assert length.shape == (1,)

            if angle is None:
                assert isinstance(inc,(float,int,np.ndarray))
                inc = np.atleast_1d(inc)
                assert inc <= 90 and inc >= -90
                sign = np.sign(inc)

                angle = (90 - np.abs(inc)) * sign
            else:
                # Assert angle between -90 and 90
                assert isinstance(angle,(float,int,np.ndarray))
                angle = np.atleast_1d(angle)
                assert angle <= 90 and angle >= -90 

            delta_z = length * np.sin(np.radians(angle))

    else:
        assert isinstance(delta_z,(float,int,np.ndarray))
        delta_z = np.atleast_1d(delta_z)
        #assert delta_z.shape == (1,)


    #Assert ge be positive
    assert isinstance(ge,(float,int,np.ndarray)) and ge>0

    #Calculate Delta P
    delta_p = 0.433 * ge * delta_z

    #Calculate P2
    p2 = p1 + delta_p

    return delta_p, p2

def kinetic_energy_change(d1=None,d2=None, ge=1,rate=None,p1=0):
    """
    Δp KE is the pressure drop resulting from a change in the velocity of the fluid between positions 1 and 2.
    It will be zero for an incompressible fluid unless the cross-sectional area of the pipe is different at the
    two positions of interest.

    Petroleum Production Systems, Economides. Chapter 7 7.2.3.2. Δp KE, the Pressure Drop Due to Kinetic Energy Change. Page 172

    """

    assert isinstance(d1,(float,int,np.ndarray)) and isinstance(d2,(float,int,np.ndarray))
    d1 = np.atleast_1d(d1)
    d2 = np.atleast_1d(d2)


    #Assert Specifi Gravity be positive
    assert isinstance(ge,(float,int,np.ndarray)) and ge>0
    ge = np.atleast_1d(ge)


    # Rate in bbl/d
    assert isinstance(rate,(float,int,np.ndarray)) and rate>0
    rate = np.atleast_1d(rate) 

    #Estimate Density in lb/ft3
    rho = 62.4 * ge

    #Estimate delta Pressure in psi
    delta_p = 1.53e-8 * np.power(rate,2) * rho * ((1/np.power(d2,4))-(1/np.power(d1,4)))

    p2 = p1 - delta_p

    return delta_p, p2

def reynolds_number(rate,rho,d,mu):
    """
    Reynolds Number where q is in bbl/d, ρ in lb m /ft 3 , D in in., and μ in cp.
    """ 
    nre = (1.48 * rate * rho) / (d * mu)

    return nre

def frictional_pressure_drop(
    rate=None, 
    epsilon=0.001,
    ge=1,
    d=None, 
    mu=1, 
    length=None):

    # Rate in bbl/d
    assert isinstance(rate,(float,int,np.ndarray)) and rate>0
    rate = np.atleast_1d(rate) 

    # pipe relative roughness
    assert isinstance(epsilon,(float,int,np.ndarray))
    epsilon = np.atleast_1d(epsilon) 

    #Assert Specifi Gravity be positive
    assert isinstance(ge,(float,int,np.ndarray)) and ge>0
    ge = np.atleast_1d(ge)

    assert isinstance(d,(float,int,np.ndarray))
    d = np.atleast_1d(d)

    assert isinstance(mu,(float,int,np.ndarray))
    mu = np.atleast_1d(mu)

    assert isinstance(length,(float,int,np.ndarray))
    length = np.atleast_1d(length)

    #Estimate Density in lb/ft3
    rho = 62.4 * ge

    #Reynolds Number
    nre = reynolds_number(rate,rho,d,mu)

    #Friction Factor
    ff = np.power((1/(-4*np.log10((epsilon/3.7065)-(5.0452/nre)*np.log10((np.power(epsilon,1.1098)/2.8257)+np.power(7.149/nre,0.8981))))),2)

    #Velocity ft/s
    u = (4*rate*5.615)/(np.pi*np.power(d/12,2)*86400)

    delta_p = (2 * ff * rho * np.power(u,2) * length)/(32.17 * (d/12) * 144)
    delta_p *= -1
    return delta_p



def incompressible_pressure_profile(
    p1=0,
    ge=1,
    epsilon=0.001,
    md=None,
    tvd=None,
    d = None,
    rate = None,
    mu=None
    ):

    assert isinstance(md,(int,float,list,np.ndarray))
    md = np.atleast_1d(md)
    assert isinstance(tvd,(int,float,list,np.ndarray))
    tvd = np.atleast_1d(tvd)
    assert isinstance(d,(int,float,list,np.ndarray))
    d = np.atleast_1d(d)
    assert isinstance(rate,(int,float))
    rate = np.atleast_1d(rate)
    assert isinstance(mu,(int,float))
    mu = np.atleast_1d(mu)
    assert isinstance(p1,(int,float))
    p1 = np.atleast_1d(p1)
    assert isinstance(ge,(int,float))
    ge = np.atleast_1d(ge)
    assert isinstance(epsilon,(int,float))
    epsilon = np.atleast_1d(epsilon)

    assert md.shape[0] == tvd.shape[0] == d.shape[0]

    n = md.shape[0]

    #Create arrays
    pressure = np.zeros(n)
    ppe = np.zeros(n)
    pke = np.zeros(n)
    pf = np.zeros(n)
    delta_p = np.zeros(n)

    pressure[0] = p1

    for i in range(1,n):

        #Potential Energy Change
        ppe[i], _ = potential_energy_change(
            z1=tvd[i-1],
            z2=tvd[i],
            ge= ge,
        )

        #Kinetic Energy Change
        pke[i], _ = kinetic_energy_change(
            d1=d[i-1],
            d2=d[i],
            rate=rate,
            ge=ge,
        )

        #Frictional Pressure drop
        pf[i] = frictional_pressure_drop(
            rate=rate, 
            epsilon=epsilon,
            ge=ge,
            d=d[i], 
            mu=mu, 
            length=np.abs(md[i-1]-md[i])
        )

        delta_p[i] = ppe[i] + pke[i] + pf[i]
        pressure[i] = pressure[i-1] + delta_p[i]

    
        # Create dataframe
    pressure_profile = pd.DataFrame({
        'md':md,
        'tvd':tvd,
        'pressure':pressure,
        'ppe': ppe,
        'pke': pke,
        'pf' : pf,
        'delta_p': delta_p
    })

    return pressure_profile




## Gas Outflow functions

def gas_pressure_profile_correlation(thp,sg,depth):
    assert isinstance(thp,(float,int,np.ndarray))
    thp = np.atleast_1d(thp)
    assert thp.ndim == 1

    assert isinstance(sg,(float,int,np.ndarray))
    sg = np.atleast_1d(sg)
    assert sg.shape == (1,)

    assert isinstance(depth,(list,float,int,np.ndarray))
    depth = np.atleast_1d(depth)
    assert sg.ndim == 1

    pwf = thp*np.exp(3.47e-5*depth)

    return pwf



def gas_pressure_profile(md, inc, thp, rate, gas_obj,di=2.99,surf_temp=80,temp_grad=1,epsilon = 0.0006, tol = 0.05, max_iter=20):
    """
    To calculate the pressure drop in a gas well, the compressibility of the fluid must be considered. When
    the fluid is compressible, the fluid density and fluid velocity vary along the pipe, and these variations
    must be included when integrating the mechanical energy balance equation.

    Petroleum Production Systems, Economides. Chapter 7 7.3. Single-Phase Flow of a Compressible, Newtonian Fluid. Page 175

    """
    # Assert the right types and shapes for input
    assert isinstance(md, (np.ndarray,pd.Series)) and md.ndim ==1
    md = np.atleast_1d(md)

    assert isinstance(inc, (int,float,np.ndarray,pd.Series))
    if isinstance(inc,np.ndarray):
        assert inc.shape == md.shape
    else:
        inc = np.full(md.shape,inc)

    angle = np.radians(90 - inc) 

    assert isinstance(thp, (int,np.int32,np.float64,float,np.ndarray)), f'{type(thp)} not accepted'
    thp = np.atleast_1d(thp)
    assert thp.shape == (1,)

    assert isinstance(gas_obj,gas) and gas_obj.pvt is not None

    assert isinstance(di, (int,float,np.ndarray))
    if isinstance(di,np.ndarray):
        assert di.shape == md.shape
    else:
        di = np.full(md.shape,di)

    assert isinstance(rate, (int,float,np.ndarray))
    rate = np.atleast_1d(rate)
    assert rate.shape == (1,)

    assert gas_obj.sg is not None

    #Create the variables

    pressure_profile = np.zeros(md.shape)
    temperature_profile = np.zeros(md.shape)
    pressure_gradient = np.zeros(md.shape)
    pressure_profile[0] = thp
    temperature_profile[0] = surf_temp

    interations = np.zeros(md.shape)

    if gas_obj.chromatography is not None:
        df_rho = gas_obj.chromatography.get_rhog(p=thp,t=surf_temp, rhog_method='real_gas')
    else:
        df_rho = gas_obj.pvt.interpolate(thp,property='rhog')

    grad_guess = df_rho['rhog'].values*(0.433/62.4)

    #Loop over depth
    for i in range(1,md.shape[0]):
        err = tol + 0.01
        dz = np.sin(angle[i])*(md[i]-md[i-1])
        gas_sg = gas_obj.sg
        it = 0
        while err>= tol and it <= max_iter:
            p_guess = grad_guess*(md[i]-md[i-1])*np.sin(angle[i]) + pressure_profile[i-1]

            #Interpolate pvt
            df_pvt = gas_obj.pvt.interpolate(p_guess)

            #Reynolds Number
            #nre = (4*28.97*gas_obj.sg*rate*14.7)/(np.pi*di[i]*df_pvt['mug'].values*10.73*520)
            nre = 20.09*(gas_sg*rate)/(di[i]*df_pvt['mug'].values)

            #Friction Factor
            friction = np.power((1/(-4*np.log10((epsilon/3.7065)-(5.0452/nre)*np.log10((np.power(epsilon,1.1098)/2.8257)+np.power(7.149/nre,0.8981))))),2)

            #Temperature
            temperature_profile[i] = dz * (temp_grad/100) + temperature_profile[i-1]

            #S
            s = (-0.0375*gas_obj.sg*dz)/(df_pvt['z'].values*(temperature_profile[i]+460))

            #Calculate next pressure by parts for easily read
            a = np.exp(-s) * np.power(pressure_profile[i-1],2)
            b = (friction*np.power(df_pvt['z'].values*(temperature_profile[i]+460)*rate,2))/(np.sin(angle[i])*np.power(di[i],5))
            c = 1 - np.exp(-s)

            p_new = np.sqrt(a - (2.685e-3*b*c))
            grad_new = (p_new - pressure_profile[i-1])/dz

            err = np.abs(grad_guess-grad_new)/grad_new
            grad_guess = grad_new
            it +=1
        
        pressure_gradient[i] = grad_new
        pressure_profile[i] = p_new
        interations[i] = it

    df_dict = {
        'pressure':pressure_profile,
        'pressure_gradient': pressure_gradient,
        'temperature': temperature_profile,
        'iterations': interations
    }

    df = pd.DataFrame(df_dict, index = md)
    pwf = pressure_profile[-1]

    return df, pwf
    """
def gas_pressure_profile_gray(md, inc, thp, rate, gas_obj,di=2.99,surf_temp=80,temp_grad=1,epsilon = 0.0006, tol = 0.05, max_iter=20):

    The Gray correlation was developed specifically for wet gas wells and is commonly used for gas wells
    producing free water and/or condensate with the gas. This correlation empirically calculates liquid
    holdup to compute the potential energy gradient and empirically calculates an effective pipe roughness
    to determine the frictional pressure gradient.

    Petroleum Production Systems, Economides. Chapter 7 7.4.3.5. The Gray Correlation page 197
    """

def gas_outflow_curve(
    md, 
    inc, 
    thp, 
    gas_obj,
    rate=None,
    min_rate=100,
    max_rate=8000,
    n_rate=20,
    di=2.99,
    surf_temp=80,
    temp_grad=1,
    epsilon = 0.0006, 
    tol = 0.05, 
    max_iter=20
    ):

    # Assert the right types and shapes for input
    assert isinstance(md, (np.ndarray,pd.Series)) and md.ndim ==1
    md = np.atleast_1d(md)

    assert isinstance(inc, (int,float,np.ndarray,pd.Series))
    if isinstance(inc,np.ndarray):
        assert inc.shape == md.shape
    else:
        inc = np.full(md.shape,inc)

    angle = np.radians(90 - inc) 

    assert isinstance(thp, (int,float,list,np.ndarray))
    thp = np.atleast_1d(thp)
    assert thp.ndim == 1

    assert isinstance(gas_obj,gas) and gas_obj.pvt is not None

    assert isinstance(di, list)

    assert isinstance(rate, (int,float,list,np.ndarray,type(None)))
    if rate is None:
        rate = np.linspace(min_rate,max_rate,n_rate)
    else:
        rate = np.atleast_1d(rate)
        assert rate.ndim == 1

    assert gas_obj.sg is not None

    pwf = np.zeros((rate.shape[0],thp.shape[0]*len(di)))
    name_list = []
    c = 0
    for p in thp:
        for d in di:
            i = 0
            for q in rate:
                _,pwf[i,c] = gas_pressure_profile(md,inc,p,q,gas_obj,surf_temp=surf_temp,temp_grad=temp_grad,di=d)
                i += 1
            c += 1
            col_name = f'thp-{p}_di-{np.mean(d)}'
            name_list.append(col_name)

    df = pd.DataFrame(pwf,columns=name_list,index=rate)

    return df

    