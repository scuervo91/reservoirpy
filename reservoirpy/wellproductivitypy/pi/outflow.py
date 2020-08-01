import numpy as np
import pandas as pd 
from ...pvtpy.black_oil import pvt, gas


## Gas Outflow function

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
            friction = np.power((1/(-4*np.log((epsilon/3.7065)-(5.0452/nre)*np.log((np.power(epsilon,1.1098)/2.8257)+np.power(7.149/nre,0.8981))))),2)

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

    