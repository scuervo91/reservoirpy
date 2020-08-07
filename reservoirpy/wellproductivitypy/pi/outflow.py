import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from ...pvtpy.black_oil import pvt, gas, oil, water

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
            assert isinstance(z1,(float,int,np.ndarray,np.int64,np.float64)) and isinstance(z2,(float,int,np.ndarray,np.int64,np.float64)), f"{type(z1)} {type(z2)}"
            z1 = np.atleast_1d(z1)
            z2 = np.atleast_1d(z2)
            #assert z1.shape == (1,) and z2.shape == (1,)
            delta_z = z1-z2

        else:
            assert isinstance(length,(float,int,np.ndarray,np.int64,np.float64)) 
            length = np.atleast_1d(length)
            #assert length.shape == (1,)

            if angle is None:
                assert isinstance(inc,(float,int,np.ndarray,np.int64,np.float64))
                inc = np.atleast_1d(inc)
                assert inc <= 90 and inc >= -90
                sign = np.sign(inc)

                angle = (90 - np.abs(inc)) * sign
            else:
                # Assert angle between -90 and 90
                assert isinstance(angle,(float,int,np.ndarray,np.int64,np.float64))
                angle = np.atleast_1d(angle)
                assert angle <= 90 and angle >= -90 

            delta_z = length * np.sin(np.radians(angle))

    else:
        assert isinstance(delta_z,(float,int,np.ndarray,np.int64,np.float64))
        delta_z = np.atleast_1d(delta_z)
        #assert delta_z.shape == (1,)


    #Assert ge be positive
    assert isinstance(ge,(float,int,np.ndarray,np.int64,np.float64)) and ge>0

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

    assert isinstance(d1,(float,int,np.ndarray,np.int64,np.float64)) and isinstance(d2,(float,int,np.ndarray,np.int64,np.float64))
    d1 = np.atleast_1d(d1)
    d2 = np.atleast_1d(d2)


    #Assert Specifi Gravity be positive
    assert isinstance(ge,(float,int,np.ndarray,np.int64,np.float64)) and ge>0
    ge = np.atleast_1d(ge)


    # Rate in bbl/d
    assert isinstance(rate,(float,int,np.ndarray,np.int64,np.float64)) and rate>0
    rate = np.atleast_1d(rate) 

    #Estimate Density in lb/ft3
    rho = 62.4 * ge

    #Estimate delta Pressure in psi
    delta_p = 1.53e-8 * np.power(rate,2) * rho * ((1/np.power(d1,4))-(1/np.power(d2,4)))

    p2 = p1 + delta_p

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
    assert isinstance(rate,(float,int,np.ndarray,np.int64,np.float64)) and rate>0
    rate = np.atleast_1d(rate) 

    # pipe relative roughness
    assert isinstance(epsilon,(float,int,np.ndarray,np.int64,np.float64))
    epsilon = np.atleast_1d(epsilon) 

    #Assert Specifi Gravity be positive
    assert isinstance(ge,(float,int,np.ndarray,np.int64,np.float64)) and ge>0
    ge = np.atleast_1d(ge)

    assert isinstance(d,(float,int,np.ndarray,np.int64,np.float64))
    d = np.atleast_1d(d)

    assert isinstance(mu,(float,int,np.ndarray,np.int64,np.float64))
    mu = np.atleast_1d(mu)

    assert isinstance(length,(float,int,np.ndarray,np.int64,np.float64))
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
        'diameter':d,
        'pressure':pressure,
        'ppe': ppe,
        'pke': pke,
        'pf' : pf,
        'delta_p': delta_p
    })

    return pressure_profile


## Gas Outflow functions

def gas_pressure_profile_correlation(thp,sg,depth):
    assert isinstance(thp,(float,int,np.ndarray,np.int64,np.float64))
    thp = np.atleast_1d(thp)
    assert thp.ndim == 1

    assert isinstance(sg,(float,int,np.ndarray,np.int64,np.float64))
    sg = np.atleast_1d(sg)
    assert sg.shape == (1,)

    assert isinstance(depth,(list,float,int,np.ndarray))
    depth = np.atleast_1d(depth)
    assert sg.ndim == 1

    pwf = thp*np.exp(3.47e-5*depth)

    return pwf



def gas_pressure_profile(
    md, 
    inc, 
    thp, 
    rate, 
    gas_obj,
    di=2.99,
    surf_temp=80,
    temp_grad=1,
    epsilon = 0.0006, 
    tol = 0.05, 
    max_iter=20):
    """
    To calculate the pressure drop in a gas well, the compressibility of the fluid must be considered. When
    the fluid is compressible, the fluid density and fluid velocity vary along the pipe, and these variations
    must be included when integrating the mechanical energy balance equation.

    Petroleum Production Systems, Economides. Chapter 7 7.3. Single-Phase Flow of a Compressible, Newtonian Fluid. Page 175

    """
    # Assert the right types and shapes for input
    assert isinstance(md, (np.ndarray,pd.Series))
    md = np.atleast_1d(md)
    assert  md.ndim ==1

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

### Multiphase Pressure Gradients

def flow_regime_plot(
    ql=None, 
    qg=None,
    d=2.99,
    sg_liquid = 1,
    surface_tension=30,
    ax=None,
    method = 'duns_ros',
    **kwargs
    ):
    """
    Plot Flow Regime from Duns and Ros Flow Regime Map
    
    Coordinates extracted from Figure7-10 Duns and Ros Flow Regime Map
    https://apps.automeris.io/wpd/

    Petroleum Production Systems, Economides. Chapter 7 7.2.3.2. Δp KE, the Pressure Drop Due to Kinetic Energy Change. Page 84

    """
    if d is not None:
        assert isinstance(d,(int,float,list,np.ndarray,pd.Series))
        d = np.atleast_1d(d)
        # Estimate Cross section Area [ft2] from diameter [in]
        a = np.power((d*0.5)/12,2)*np.pi

    if ql is not None:
        assert isinstance(ql,(int,float,list,np.ndarray,pd.Series))
        ql = np.atleast_1d(ql)

        #Liquid velocity. Convert bbl/d to ft3/s then divide area. Result velocity in ft/s
        usl = (ql * 5.616 * (1/86400))/a
        #Calculate the dimensionless numbers for each phase
        nvl = 1.938 * usl * np.power((sg_liquid*62.4)/surface_tension,0.25)

    if qg is not None:
        assert isinstance(ql,(int,float,list,np.ndarray,pd.Series))
        qg = np.atleast_1d(qg)

        #Gas velocity. Convert ft3/d to ft3/s then divide area. Result velocity in ft/s
        usg = (qg * (1/86400))/a
        nvg = 1.938 * usg * np.power((sg_liquid*62.4)/surface_tension,0.25)

    if method == 'duns_ros':
        fax= ax or plt.gca()
        region_1_2 = np.array([
            [1.1753722651306362, 0.1082636733874053],
            [1.1913061720030635, 0.16102620275609392],
            [1.3268047497147244, 0.23950266199874834],
            [1.4777148689707504, 0.35154183187529914],
            [1.7604108438655526, 0.5228664844415476],
            [2.1544346900318843, 0.7880462815669913],
            [2.8585141796844757, 1.2358165955824107],
            [3.545745842465605, 1.790084628235539],
            [5.529553425383406, 3.2470894518548166],
            [8.507942799627454, 5.512889788770675],
            [16.68100537200059, 11.566937549363251],
            [29.76351441631322, 20.43359717856943],
            [61.58482110660267, 39.079952122756026],
            [41.11829402435837, 27.703123342457815],
            [79.53985507023424, 48.93900918477497],
        ])

        region_2_t = np.array([
            [53.10631887314356, 0.10543589908346815],
            [59.146605445917515, 0.18139306939110614],
            [66.7669293918757, 0.36097012876068046],
            [80.61813527211957, 0.7674630429274295],
            [104.12232560483065, 1.5475873545578884],
            [141.92103954525945, 2.7338936055226313],
            [270.8622850933671, 5.9684569951223105],
            [204.14630347954724, 4.230939172613499],
            [340.53655850163904, 7.674630429274299],
            [503.2159359259993, 12.195704601594414],
            [714.1692874235849, 18.380944176677932],
            [922.3851039358485, 23.324701361610806],
        ])

        region_t_3 = np.array([
            [92.23851039358486, 0.10684043121253317],
            [97.34285811778867, 0.15475873545578891],
            [105.53385749880759, 0.24269312356542563],
            [115.96514767613999, 0.41204298882016666],
            [136.30221830031346, 0.7278953843983147],
            [183.29807108324394, 1.2358165955824107],
            [263.6650898730361, 2.271547585601246],
            [364.25331154496416, 4.120429888201667],
            [531.0631887314356, 6.995642156712631],
            [714.1692874235849, 11.264816923358868],
            [947.5632026539927, 18.139306939110632],
        ])

        fax.plot(region_1_2[:,0],region_1_2[:,1], color='black',linestyle='--')
        fax.plot(region_2_t[:,0],region_2_t[:,1], color='black',linestyle='--')
        fax.plot(region_t_3[:,0],region_t_3[:,1], color='black',linestyle='--')
        fax.set_ylabel('Nvl')
        fax.set_ylabel('Nvg')
        fax.set_title('Duns and Ros Flow Regime Map')
        fax.set_xlim([0.1,1000])
        fax.set_ylim([0.1,100])
        annot = kwargs.pop('ann',True)
        font = kwargs.pop('fontsize',8)

        if annot:
            fax.annotate(
                f"Region I \n Bubble Flow or \n low-velocity slug flow",
                xy = (0.2,0.15),
                xycoords='data',
                xytext=(0, 0), 
                textcoords='offset points',
                bbox={'boxstyle':'round', 'fc':'0.8'},
                fontsize = font
            )

            fax.annotate(
                f"Region II \n High-velocity Flow or \n churn flow",
                xy = (2,0.15),
                xycoords='data',
                xytext=(0, 0), 
                textcoords='offset points',
                bbox={'boxstyle':'round', 'fc':'0.8'},
                fontsize = font
            )

            fax.annotate(
                f"Region III \n Annular Flow Pattern",
                xy = (300,0.15),
                xycoords='data',
                xytext=(0, 0), 
                textcoords='offset points',
                bbox={'boxstyle':'round', 'fc':'0.8'},
                fontsize = font
            )
        
        if ql is not None and qg is not None:
            fax.scatter(nvg,nvl,color='blue',marker = "^")


    if method == 'taitel_dukler':
        fax= ax or plt.gca()

        region_E = np.array([
            [14.977474763452001, 0.0022033318988979545],
            [14.977474763452001, 0.006595844345274293],
            [14.977474763452001, 0.04746934676639568],
            [14.777148689707504, 0.9165263295637442],
            [14.977474763452001, 6.87270243904312],
            [14.977474763452001, 15.857064005032758]
        ])

        region_A = np.array([
            [0.08858667904100832, 0.0022372323125884317],
            [0.08858667904100832, 0.005091596044287256],
            [0.0986624843178949, 0.018460289732281962],
            [0.11137395078578621, 0.04142593768347061],
            [0.1326804749714725, 0.08679099331751502],
            [0.1668100537200059, 0.18431459769950134],
            [0.21256187881919958, 0.3275265038954424],
            [0.30575961084169306, 0.695276382058884],
            [0.46415888336127775, 1.2691784682206282],
            [0.7336637748600019, 2.019816384578137],
            [0.9223851039358476, 2.412109197346714]
        ])
        region_B = np.array([
            [0.028585141796844758, 3.4805610999729812],
            [0.0531063188731435, 3.5220947122633963],
            [0.08623280529014943, 3.517016970779084],
            [0.24649769667586238, 3.2292570594299215],
            [0.8978760230238888, 2.4455928433916867],
            [2.0971883035581533, 1.7556200043179786],
            [5.239601353002639, 4.20919831000811],
            [10.412232560483055, 7.572933314656229],
            [14.579502008614657, 10.657087726496014],
        ])
        region_D = np.array([
            [0.26366508987303583, 0.44861391200434203],
            [0.30575961084169306, 0.4018483957905594],
            [0.4398198780581129, 0.2288467215238852],
            [0.5032159359259996, 0.16920697751727592],
            [0.5835551032264551, 0.11058672774921392],
            [0.6676692939187563, 0.05647578739286295],
            [0.6951927961775606, 0.03743162248826758],
            [0.7536903980898542, 0.02284801683862376],
            [0.7639077845044221, 0.015565548854263186],
            [0.7436096708208817, 0.011357807043115235],
            [0.7847599703514607, 0.006933286608265855],
            [0.7536903980898542, 0.0027304200384003397],
            [0.7436096708208817, 0.002162999360197944],
        ])

        fax.plot(region_A[:,0],region_A[:,1], color='black',linestyle='--')
        fax.plot(region_B[:,0],region_B[:,1], color='black',linestyle='--')
        fax.plot(region_D[:,0],region_D[:,1], color='black',linestyle='--')
        fax.plot(region_E[:,0],region_E[:,1], color='black',linestyle='--')
        fax.set_ylabel('Usg [m/s]')
        fax.set_ylabel('Usl [m/s]')
        fax.set_title('Taitel-Dukler flow regime map')
        fax.set_xlim([0.01,100])
        fax.set_ylim([0.001,10])
        if ql is not None and qg is not None:
            fax.scatter(usg*0.3048,usl*0.3048,color='blue',marker = "^")


    fax.set_yscale('log')
    fax.set_xscale('log')

def hb_correlation(
    p,  #Pressure [psi]
    t, #Temperature [F]
    ql, # Liquid Flow [bbl/d]
    qg, # gas flow [kscfd]
    ten_liquid, #Surface tension dyne/cm2
    rho_liquid, # density lb/ft3
    mu_liquid, # Viscosity [cp]
    z, # Gas compressibility Factor
    di # Diameter
):
    """
    The modified Hagedorn and Brown method (mH-B) is an empirical two-phase flow correlation based
    on the original work of Hagedorn and Brown (1965). The heart of the Hagedorn-Brown method is a
    correlation for liquid holdup; the modifications of the original method include using the no-slip holdup
    when the original empirical correlation predicts a liquid holdup value less than the no-slip holdup and
    the use of the Griffith correlation (Griffith and Wallis, 1961) for the bubble flow regime.

    Petroleum Production Systems, Economides. Chapter 7 7.4.3.1. The Modified Hagedorn and Brown Method  Page 187

    """
    usl = (ql * 5.615)/(0.0278 * 86400)
    usg = (4*qg*1000*z*(460+t*14.7)/(p*520*np.pi*np.power(di/12,2)) 
    
    #Mixure Velocity
    um = usl + usg 
    



def two_phase_pressure_profile(
    depth,
    thp,
    oil_rate,
    gas_rate,
    bsw,
    oil_obj,
    gas_obj,
    water_obj, 
    epsilon=0.0006, 
    surface_temperature=80, 
    temperature_gradient=1,  
    di=2.99, 
    tol=0.05,
    max_iter = 20,
    method = 'hagedorn_brown'
):



    # Assert the right types and shapes for input
    assert isinstance(depth, (np.ndarray,pd.Series,list))
    depth = np.atleast_1d(depth)
    assert depth.ndim == 1

    assert isinstance(thp, (int,np.int32,np.float64,float,np.ndarray)), f'{type(thp)} not accepted'
    thp = np.atleast_1d(thp)
    assert thp.shape == (1,)

    assert isinstance(oil_rate, (int,np.int32,np.float64,float,np.ndarray)), f'{type(thp)} not accepted'
    oil_rate = np.atleast_1d(oil_rate)
    assert oil_rate.shape == (1,)

    assert isinstance(gas_rate, (int,np.int32,np.float64,float,np.ndarray)), f'{type(thp)} not accepted'
    gas_rate = np.atleast_1d(gas_rate)
    assert gas_rate.shape == (1,)

    assert isinstance(bsw, (int,np.int32,np.float64,float,np.ndarray)), f'{type(thp)} not accepted'
    bsw = np.atleast_1d(bsw)
    assert bsw.shape == (1,)

    assert isinstance(gas_obj,gas) and gas_obj.pvt is not None
    assert isinstance(oil_obj,oil) and oil_obj.pvt is not None
    assert isinstance(water_obj,water) and water_obj.pvt is not None

    if isinstance(di,(np.ndarray,pd.Series,list)):
        assert di.shape == depth.shape
    elif isinstance(di,(int,float)):
        di = np.full(depth.shape,di)

    assert isinstance(epsilon,(int,float,np.ndarray))
    epsilon = np.atleast_1d(epsilon)
    assert epsilon.shape == (1,)

    assert isinstance(surface_temperature,(int,float,np.ndarray))
    surface_temperature = np.atleast_1d(surface_temperature)

    assert isinstance(temperature_gradient,(int,float,np.ndarray))
    temperature_gradient = np.atleast_1d(temperature_gradient)

    #Start
    qo = oil_rate * (5.615/86400)
    ql = qo / (1-bsw)
    qw = ql * bsw

    pressure_profile = np.zeros(depth.shape)
    pressure_profile[0] = thp
    pressure_gradient = np.zeros(depth.shape)
    iterations = np.zeros(depth.shape)

    temperature_profile = np.abs(depth[0] - depth) * temperature_gradient + surface_temperature

    #Initials Densities
    rho_oil_i = oil_obj.pvt.interpolate(thp,property = 'rho').iloc[0,0]
    rho_water_i = water_obj.pvt.interpolate(thp,property = 'rho').iloc[0,0]
    rho_l = rho_oil_i * (1-bsw) + rho_oil_i * bsw 

    pressure_gradient[0] = rho_l * (0.433/62.4)
    for i in range(1,depth.shape[0]):
        err = tol + 0.01
        it = 0
        grad_guess = pressure_gradient[i-1]

        while err>= tol and it <= max_iter:
            p_guess = grad_guess * np.abs(depth[i] - depth[i-1]) + pressure_profile[i-1]
            
            #Interpolate pvt
            gas_pvt_guess = gas_obj.pvt.interpolate(p_guess)
            oil_pvt_guess = oil_obj.pvt.interpolate(p_guess)
            water_pvt_guess = water_obj.pvt.interpolate(p_guess)

            if method == 'hagedorn_brown':
                grad_new = hb_correlation()
            elif method == 'beggs_brill':
                grad_new = bb_correlation()
            elif method == 'gray':
                grad_new = bb_correlation()
            else:
                raise ValueError('No method has been chosen')

            err =  abs(grad_guess-grad_new)/grad_new
            grad_guess = grad_new
            it += 1

        pressure_gradient[i] = grad_new 
        pressure_profile[i] = p_guess
        iterations[i] = it

    df_dict = {
        'pressure':pressure_profile,
        'pressure_gradient': pressure_gradient,
        'temperature': temperature_profile,
        'iterations': iterations
    }

    df = pd.DataFrame(df_dict, index = depth)
    pwf = pressure_profile[-1]

    return df, pwf

    

    


