import numpy as np 
import pandas as pd 
from scipy.interpolate import interp1d
from .petroequations import vshale_gr, vshale_dn, phi_rho, phie, sw, perm, flow_capacity, phia, sw_pnn


def petrophysics(logs,dfrom,dto,
                vshale_gr_kw=None,       
                vshale_dn_kw=None,
                phi_rho_kw=None,    #[DenM,DenF,DenColumn,RhoPhiColumn]
                phie_kw=None,      #[RhophiColumn,NtrColumn,VshColumn,[ListMethod],[listNames]]
                sw_kw=None,        #[a,m,n,Rw,phicolumn,rtcolumn,Vshcolumn,Rsh=4.0,alpha=0.3,[ListMethod],[listNames]]
                perm_kw=None,      #[phiecolumn,swcolum,autor,fluid,[ListMethod],[listNames]]
                flag_kw=None,   #[phicolumn,phicut,vshcol,vshcut,swcol,swcut,kcol,kcut,paycol]
                kh_kw=None,       #[h,kcol,paycol,khcol]
                sw_pnn_kw = None, 
                return_partial = False
                ):
    """petrophysics [summary]

    Parameters
    ----------
    logs : [type]
        [description]
    dfrom : [type]
        [description]
    dto : [type]
        [description]
    vshale_gr_kw : [type], optional
        [description], by default None
    vshale_dn_kw : [type], optional
        [description], by default None
    phi_rho_kw : [type], optional
        [description], by default None
    return_partial : bool, optional
        [description], by default False

    Returns
    -------
    [type]
        [description]
    """
    logf=logs[(logs.index >= dfrom) & (logs.index<=dto)].copy()
    new_cols = []

    if sw_pnn_kw is not None:
        #name for the new columns
        sw_col = sw_pnn_kw.pop('sw_col_name','sw_pnn')

        #Name for input columns
        phie_name = sw_pnn_kw.pop('phie_name','phie')
        vsh_name = sw_pnn_kw.pop('vsh_name','vsh')
        sigma_name = sw_pnn_kw.pop('sigma_name','sigma')
        sighy = sw_pnn_kw.pop('sighy',20)
        sigsh = sw_pnn_kw.pop('sigsh',35)
        sigmam = sw_pnn_kw.pop('sigmam',None)
        if isinstance(sigmam,str):
            _sigmam = logf[sigmam]
        elif isinstance(sigmam,(int,float,type(None))):
            _sigmam = sigmam
        sigw = sw_pnn_kw.pop('sigw',None)

        ws = sw_pnn_kw.pop('ws',None)

        logf[sw_col] = sw_pnn(logf[phie_name],logf[vsh_name], logf[sigma_name],sighy,sigsh,sigw=sigw, sigmam = _sigmam,ws=ws)
        new_cols.append(sw_col)
    if vshale_gr_kw is not None:
        #name for the new columns
        vsh_col_name = vshale_gr_kw.pop('vsh_col_name','vsh_gr')
        vsh_type = vshale_gr_kw.pop('type','linear')
        new_cols.append(vsh_col_name)
        
        #Name for input columns
        gr_col_name = vshale_gr_kw.pop('gr_name',None)
        gr_sand = vshale_gr_kw.pop('gr_sand',None)
        gr_shale = vshale_gr_kw.pop('gr_shale',None)
        logf[vsh_col_name]=vshale_gr(logf[gr_col_name],gr_sand,gr_shale,type=vsh_type)
    
    if vshale_dn_kw is not None:
        #name for the new columns
        vsh_col_name = vshale_dn_kw.pop('vsh_col_name','vsh_dn')
        new_cols.append(vsh_col_name)
        
        #Name for input columns
        rho_col_name = vshale_dn_kw.pop('rho_name',None)
        ntr_col_name = vshale_dn_kw.pop('ntr_name',None)
        logf[vsh_col_name] = vshale_dn(logf[rho_col_name], logf[ntr_col_name], **vshale_dn_kw)
        
    if phi_rho_kw is not None:
        #name for the new columns
        phi_rho_col_name = phi_rho_kw.pop('phi_rho_name','rho_phi')
        new_cols.append(phi_rho_col_name)
        
        #Name for input columns
        rho_col_name = phi_rho_kw.pop('rho_name',None)
        logf[phi_rho_col_name]=phi_rho(logf[rho_col_name], **phi_rho_kw)
    
    if phie_kw is not None:
        #name for the new columns
        phie_avg_col_name = phie_kw.pop('phie_avg_col_name','phie_avg')
        phie_rho_col_name = phie_kw.pop('phie_rho_col_name','phie_rho')
        phie_ntr_col_name = phie_kw.pop('phie_ntr_col_name','phie_ntr')
        phi_avg_col_name = phie_kw.pop('phi_avg_col_name','phia')
        
        #Name for input columns
        phi_rho_col_name = phie_kw.pop('phi_rho_name',None)
        ntr_col_name = phie_kw.pop('ntr_name',None)
        vsh_col_name = phie_kw.pop('vsh_name',None)
        if (phi_rho_col_name is not None) & (ntr_col_name is not None):
            logf[phi_avg_col_name] = phia(logf[phi_rho_col_name],logf[ntr_col_name], **phie_kw)
            logf[phie_avg_col_name]=phie(logf[phi_avg_col_name],logf[vsh_col_name])
            logf[phie_rho_col_name]=phie(logf[phi_rho_col_name],logf[vsh_col_name])
            logf[phie_ntr_col_name]=phie(logf[ntr_col_name],logf[vsh_col_name])
            new_cols.extend([phi_avg_col_name,phie_avg_col_name,phie_rho_col_name,phie_ntr_col_name])
        elif phi_rho_col_name is not None:
            logf[phie_rho_col_name]=phie(logf[phi_rho_col_name],logf[vsh_col_name])
            new_cols.append(phie_rho_col_name)
        elif ntr_col_name is not None:
            logf[phie_ntr_col_name]=phie(logf[ntr_col_name],logf[vsh_col_name])
            new_cols.append(phie_ntr_col_name)
    
    if sw_kw is not None:
        #Name for input columns
        rt_col_name = sw_kw.pop('rt_name',None)
        phi_col_name = sw_kw.pop('phi_name',None)
        vsh_col_name = sw_kw.pop('vsh_name',None)
        sw_kw['vsh_curve'] = logf[vsh_col_name] if vsh_col_name!=None else None
        rw = sw_kw.pop('rw', None)
        methods = sw_kw.pop('methods',['archie'])
        
        #name for the new columns. List of curve names depending on the method
        sw_cols_name = sw_kw.pop('sw_cols_name',methods) 
        
        for i,method in enumerate(methods):
            logf['sw_' + sw_cols_name[i]] = sw(logf[rt_col_name],logf[phi_col_name],rw,method=method, **sw_kw)
            new_cols.append(f'sw_{sw_cols_name[i]}')

    if perm_kw is not None:
        #Name for input columns
        phi_col_name = perm_kw.pop('phi_name',None)
        swir = perm_kw.pop('swir', None)
        authors = perm_kw.pop('authors',['timur'])
        fluid = perm_kw.pop('fluid','oil')

        
        #name for the new columns. List of curve names depending on the method
        perm_cols_name = sw_kw.pop('perm_cols_name',authors) 
        
        for i,author in enumerate(authors):
            logf['k_' + perm_cols_name[i]] = perm(logf[phi_col_name],swir,author=author,fluid=fluid)
            new_cols.append(perm_cols_name[i])

    
    if flag_kw is not None:
        #name for the new columns.
        sand_flag_col_name = flag_kw.pop('sand_flag_name','sand_flag')
        reservoir_flag_col_name = flag_kw.pop('reservoir_flag_name','reservoir_flag')
        pay_flag_col_name = flag_kw.pop('pay_flag_name','pay_flag')
        
        #Name for input columns
        vsh_col_name = flag_kw.pop('vsh_name',None)
        phi_col_name = flag_kw.pop('phi_name',None)
        sw_col_name = flag_kw.pop('sw_name',None)
        
        #Name for cutoffs
        vsh_cutoff = flag_kw.pop('vsh_cutoff',0)
        phi_cutoff = flag_kw.pop('phi_cutoff',0)
        sw_cutoff = flag_kw.pop('sw_cutoff',0)
        
        #whichs flags
        method = flag_kw.pop('which',None)
        
        if method=='pay':
            logf[sand_flag_col_name] = (logf[vsh_col_name] <= vsh_cutoff)*1
            logf[reservoir_flag_col_name] = (logf[phi_col_name] >= phi_cutoff)*logf[sand_flag_col_name]
            logf[pay_flag_col_name] = (logf[sw_col_name] <= sw_cutoff)*logf[reservoir_flag_col_name]
            new_cols.extend([sand_flag_col_name,reservoir_flag_col_name,pay_flag_col_name])
        elif method=='reservoir':
            logf[sand_flag_col_name] = (logf[vsh_col_name] <= vsh_cutoff)*1
            logf[reservoir_flag_col_name] = (logf[phi_col_name] >= phi_cutoff)*logf[sand_flag_col_name]
            new_cols.extend([sand_flag_col_name,reservoir_flag_col_name])
        elif method=='sand':
            logf[sand_flag_col_name] = (logf[vsh_col_name] <= vsh_cutoff)*1
            new_cols.append(sand_flag_col_name)
            
    if kh_kw is not None:
        #name for the new columns.
        kh_col_name = kh_kw.pop('kh_name','kh')
        kh_norm_col_name = kh_kw.pop('khnorm_name','kh_norm')
        
        #Name for input columns
        perm_col_name = kh_kw.pop('perm_name',None)
        pay_col_name = kh_kw.pop('pay_name','pay_flag')
        
        #
        h=np.mean(np.diff(logf.index))
        logf[kh_col_name],logf[kh_norm_col_name]=flow_capacity(h,logf[perm_col_name],logf[pay_col_name])
        new_cols.extend([kh_norm_col_name,kh_col_name])
                                                                                                                                                                                                   
    if return_partial:
        return logf
    else:
        log_merged = logs.merge(logf[new_cols], how='left', left_index=True, right_index=True)
        return log_merged