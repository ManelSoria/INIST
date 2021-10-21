function [ret] = INIST(varargin)
% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2021
%
% Units: T(K), p(bar), h and u: kJ/kg, v: m^3/kg, rho: kg/m^3 s: kJ/kgK,
% a: m/s, cv and cp: kJ/kgK, JT: bar/K, mu: Pa.s, k: W/mK, MM: kg/mol
% SF: N.m
% 1st argument: substance name
%               'Database' to return the list of database elements
% 2nd and remaining arguments: 
%  critical temperature     'tcrit' 
%  critical pressure        'pcrit'
%  critical volume          'vcrit' 
%  molecular mass           'MM'
%  saturation temperature   'tsat_p', p  
%  saturation pressure      'psat_t', T 
%  saturated liquid properties as a function of pressure
%     volume          'vl_p' , p 
%     energy          'ul_p' , p
%     enthalpy        'hl_p' , p
%     entropy         'sl_p' , p
%     specific heat coeff at constant volume:
%                     'cvl_p', p
%     specific heat coeff at constant pressure:
%                     'cpl_p', p
%     sound speed     'al_p' , p
%     viscosity       'mul_p', p
%     density         'rl_p' , p
%     conductivity    'kl_p' , p
%  saturated vapour properties as a function of pressure
%     volume          'vv_p' , p 
%     energy          'uv_p' , p
%     enthalpy        'hv_p' , p
%     entropy         'sv_p' , p
%     specific heat coeff at constant volume:
%                     'cvv_p', p
%     specific heat coeff at constant pressure:
%                     'cpv_p', p
%     sound speed     'av_p' , p
%     viscosity       'muv_p', p
%     density         'rv_p' , p
%     conductivity    'kv_p' , p
%  saturated liquid properties as a function of temperature
%     volume          'vl_t' , t 
%     energy          'ul_t' , t
%     enthalpy        'hl_t' , t
%     entropy         'sl_t' , t
%     specific heat coeff at constant volume:
%                     'cvl_t', t
%     specific heat coeff at constant pressure:
%                     'cpl_t', t
%     sound speed     'al_t' , t
%     viscosity       'mul_t', t
%     density         'rl_t' , t
%     conductivity    'kl_t' , t
%  saturated vapour properties as a function of temperature
%     volume          'vv_t' , t 
%     energy          'uv_t' , t
%     enthalpy        'hv_t' , t
%     entropy         'sv_t' , t
%     specific heat coeff at constant volume:
%                     'cvv_t', t
%     specific heat coeff at constant pressure:
%                     'cpv_t', t
%     sound speed     'av_t' , t
%     viscosity       'muv_t', t
%     density         'rv_t' , t
%     conductivity    'kv_t' , t
%  saturated vapour-liquid properties as a function of temperature and
%       quality
%     volume          'v_tx' , t, x 
%     energy          'u_tx' , t, x 
%     enthalpy        'h_tx' , t, x 
%     entropy         's_tx' , t, x 
%     specific heat coeff at constant volume:
%                     'cv_tx' , t, x 
%     specific heat coeff at constant pressure:
%                     'cp_tx' , t, x 
%     sound speed     'a_tx' , t, x 
%     viscosity       'mu_tx' , t, x 
%     density         'r_tx' , t, x 
%     conductivity    'k_tx' , t, x 
%  saturated vapour-liquid properties as a function of pressure and
%       quality
%     volume          'v_px' , p, x 
%     energy          'u_px' , p, x 
%     enthalpy        'h_px' , p, x 
%     entropy         's_px' , p, x 
%     specific heat coeff at constant volume:
%                     'cv_px' , p, x 
%     specific heat coeff at constant pressure:
%                     'cp_px' , p, x 
%     sound speed     'a_px' , p, x 
%     viscosity       'mu_px' , p, x 
%     density         'r_px' , p, x 
%     conductivity    'k_px' , p, x 
%  non-saturated properties as a function of pressure and temperature
%     volume          'v_pt' , p , t
%     energy          'u_pt' , p , t
%     enthaply        'h_pt' , p , t
%     entrophy        's_pt' , p , t
%     specific heat coeff at constant volume:
%                     'cv_pt', p , t
%     specific heat coeff at constant pressure:
%                     'cp_pt', p , t
%     sound speed     'a_pt' , p , t
%     viscosity       'mu_pt', p , t
%     density         'r_pt' , p , t
%     conductivity    'k_pt',  p , t
%  temperature as a function of ...
%     pressure and entropy 't_ps', p ,s  
%     pressure and enthalpy 't_ph', p ,h  
%     
%  special functions:
%       'minp'        returns the minimum isobar available
%       'maxp'        idem max isobar
%       'mint'        idem minimum temperature 
%       'maxt'        idem maximum temperature 
%       'isobars'     returns a vector with the available isobars
% 
%       'add_p' , p, filename     Disabled
%                                 downloads and stores a new isobar.
%                                 filename is optional. If given, it is 
%                                 the name of the file where the IND
%                                 data is stored needs to be
%                                 specified, so that the isobar
%                                 data can be updated. This will add
%                                 the new isobar to the database file.
% 

global IND

if 7~=exist('Database','dir')
    path = fileparts(which(mfilename));
    addpath(genpath(path));
end

if strcmp(varargin{1},'Database') % Return a list of species in the database 
    databasepath = [path '\Database\*.mat'];
    % Adapts path to the OS
    databasepath = osi(databasepath);
    info  = dir(databasepath);
    ret = {info.name};
    ret = strrep(ret,'.mat','');
    return
end

try % load the species needed
    if isempty(IND) || ~isfield(IND,varargin{1})  
        set = load(varargin{1});
        IND.(varargin{1}) = set.(varargin{1});
    end
catch
    error('%s not found',varargin{1})
end


dat = IND.(varargin{1});
prop = lower(varargin{2}); 

switch prop
    case 'mm'
        ret=dat.MM;
    case 'tcrit'
        ret=dat.Tsat(end);
    case 'pcrit'
        ret=dat.Psat(end);
    case 'vcrit'
        ret=dat.vl(end);
    case 'isobars'
        ret=[];
        for i=1:length(dat.isoP)
            ret(i)=dat.isoP{i}.P;
        end
    case 'tsat_p' 
        p=varargin{3};
        check_satp(p);
        ret = interp1(dat.Psat,dat.Tsat,p);
    case 'psat_t'
        T=varargin{3};
        check_satT(T);
        ret = interp1(dat.Tsat,dat.Psat,T);   
        
%  saturated liquid properties as a function of pressure    
    case 'vl_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.vl,p);        
    case 'ul_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.ul,p);
    case 'hl_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.hl,p);
    case 'sl_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.sl,p);
    case 'cvl_p', p=varargin{3}; check_undefined_p(p); check_satp(p);ret=interp1(dat.Psat,dat.cvl,p);
    case 'cpl_p', p=varargin{3}; check_undefined_p(p); check_satp(p);ret=interp1(dat.Psat,dat.cpl,p);
    case 'al_p',  p=varargin{3}; check_undefined_p(p); check_satp(p);ret=interp1(dat.Psat,dat.al,p);
    case 'mul_p', p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.mul,p);
    case 'rl_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.rl,p);
    case 'kl_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.kl,p);
        
%  saturated vapour properties as a function of pressure        
    case 'vv_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.vv,p);        
    case 'uv_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.uv,p);
    case 'hv_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.hv,p);
    case 'sv_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.sv,p);
    case 'cvv_p', p=varargin{3}; check_undefined_p(p); check_satp(p);ret=interp1(dat.Psat,dat.cvv,p);
    case 'cpv_p', p=varargin{3}; check_undefined_p(p); check_satp(p);ret=interp1(dat.Psat,dat.cpv,p);
    case 'av_p',  p=varargin{3}; check_undefined_p(p); check_satp(p);ret=interp1(dat.Psat,dat.av,p);    
    case 'muv_p', p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.muv,p);
    case 'rv_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.rv,p);
    case 'kv_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.kv,p);    
    
%  saturated liquid properties as a function of temperature        
    case 'vl_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.vl,t);        
    case 'ul_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.ul,t);
    case 'hl_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.hl,t);
    case 'sl_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.sl,t);
    case 'cvl_t', t=varargin{3}; check_undefined_t(t); check_satT(t);ret=interp1(dat.Tsat,dat.cvl,t);
    case 'cpl_t', t=varargin{3}; check_undefined_t(t); check_satT(t);ret=interp1(dat.Tsat,dat.cpl,t);
    case 'al_t',  t=varargin{3}; check_undefined_t(t); check_satT(t);ret=interp1(dat.Tsat,dat.al,t);    
    case 'mul_t', t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.mul,t);
    case 'rl_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.rl,t);
    case 'kl_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.kl,t);   

        
%  saturated vapour properties as a function of temperature
    case 'vv_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.vv,t);        
    case 'uv_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.uv,t);
    case 'hv_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.hv,t);
    case 'sv_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.sv,t);
    case 'cvv_t', t=varargin{3}; check_undefined_t(t); check_satT(t);ret=interp1(dat.Tsat,dat.cvv,t);
    case 'cpv_t', t=varargin{3}; check_undefined_t(t); check_satT(t);ret=interp1(dat.Tsat,dat.cpv,t);
    case 'av_t',  t=varargin{3}; check_undefined_t(t); check_satT(t);ret=interp1(dat.Tsat,dat.av,t);
    case 'muv_t', t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.muv,t);
    case 'rv_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.rv,t);
    case 'kv_t',  t=varargin{3}; check_satT(t);ret=interp1(dat.Tsat,dat.kv,t);  
        
%  saturated  properties as a function of temperature and quality
    case 'v_tx',  t=varargin{3};x=varargin{4}; check_satT(t);vap=interp1(dat.Tsat,dat.vv,t);liq=interp1(dat.Tsat,dat.vl,t);ret = x*(vap-liq)+liq;     
    case 'u_tx',  t=varargin{3};x=varargin{4}; check_satT(t);vap=interp1(dat.Tsat,dat.uv,t);liq=interp1(dat.Tsat,dat.ul,t);ret = x*(vap-liq)+liq;
    case 'h_tx',  t=varargin{3};x=varargin{4}; check_satT(t);vap=interp1(dat.Tsat,dat.hv,t);liq=interp1(dat.Tsat,dat.hl,t);ret = x*(vap-liq)+liq;
    case 's_tx',  t=varargin{3};x=varargin{4}; check_satT(t);vap=interp1(dat.Tsat,dat.sv,t);liq=interp1(dat.Tsat,dat.sl,t);ret = x*(vap-liq)+liq;
    case 'cv_tx', t=varargin{3};x=varargin{4}; check_satT(t);check_undefined_t(t);vap=interp1(dat.Tsat,dat.cvv,t);liq=interp1(dat.Tsat,dat.cvl,t);ret = x*(vap-liq)+liq;
    case 'cp_tx', t=varargin{3};x=varargin{4}; check_satT(t);check_undefined_t(t);vap=interp1(dat.Tsat,dat.cpv,t);liq=interp1(dat.Tsat,dat.cpl,t);ret = x*(vap-liq)+liq;
    case 'a_tx',  t=varargin{3};x=varargin{4}; check_satT(t);check_undefined_t(t);vap=interp1(dat.Tsat,dat.av,t);liq=interp1(dat.Tsat,dat.al,t);ret = x*(vap-liq)+liq;
    case 'mu_tx', t=varargin{3};x=varargin{4}; check_satT(t);vap=interp1(dat.Tsat,dat.muv,t);liq=interp1(dat.Tsat,dat.mul,t);ret = x*(vap-liq)+liq;
    case 'r_tx',  t=varargin{3};x=varargin{4}; check_satT(t);vap=interp1(dat.Tsat,dat.rv,t);liq=interp1(dat.Tsat,dat.rl,t);ret = x*(vap-liq)+liq;
    case 'k_tx',  t=varargin{3};x=varargin{4}; check_satT(t);vap=interp1(dat.Tsat,dat.kv,t);liq=interp1(dat.Tsat,dat.kl,t);ret = x*(vap-liq)+liq; 
        
%   saturated  properties as a function of pressure and quality
    case 'v_px',  p=varargin{3};x=varargin{4}; check_satp(p);vap=interp1(dat.Tsat,dat.vv,p);liq=interp1(dat.Tsat,dat.vl,p);ret = x*(vap-liq)+liq;     
    case 'u_px',  p=varargin{3};x=varargin{4}; check_satp(p);vap=interp1(dat.Tsat,dat.uv,p);liq=interp1(dat.Tsat,dat.ul,p);ret = x*(vap-liq)+liq;
    case 'h_px',  p=varargin{3};x=varargin{4}; check_satp(p);vap=interp1(dat.Tsat,dat.hv,p);liq=interp1(dat.Tsat,dat.hl,p);ret = x*(vap-liq)+liq;
    case 's_px',  p=varargin{3};x=varargin{4}; check_satp(p);vap=interp1(dat.Tsat,dat.sv,p);liq=interp1(dat.Tsat,dat.sl,p);ret = x*(vap-liq)+liq;
    case 'cv_px', p=varargin{3};x=varargin{4}; check_satp(p);check_undefined_p(p);vap=interp1(dat.Tsat,dat.cvv,p);liq=interp1(dat.Tsat,dat.cvl,p);ret = x*(vap-liq)+liq;
    case 'cp_px', p=varargin{3};x=varargin{4}; check_satp(p);check_undefined_p(p);vap=interp1(dat.Tsat,dat.cpv,p);liq=interp1(dat.Tsat,dat.cpl,p);ret = x*(vap-liq)+liq;
    case 'a_px',  p=varargin{3};x=varargin{4}; check_satp(p);check_undefined_p(p);vap=interp1(dat.Tsat,dat.av,p);liq=interp1(dat.Tsat,dat.al,p);ret = x*(vap-liq)+liq;
    case 'mu_px', p=varargin{3};x=varargin{4}; check_satp(p);vap=interp1(dat.Tsat,dat.muv,p);liq=interp1(dat.Tsat,dat.mul,p);ret = x*(vap-liq)+liq;
    case 'r_px',  p=varargin{3};x=varargin{4}; check_satp(p);vap=interp1(dat.Tsat,dat.rv,p);liq=interp1(dat.Tsat,dat.rl,p);ret = x*(vap-liq)+liq;
    case 'k_px',  p=varargin{3};x=varargin{4}; check_satp(p);vap=interp1(dat.Tsat,dat.kv,p);liq=interp1(dat.Tsat,dat.kl,p);ret = x*(vap-liq)+liq;

        
    case {'v_pt','u_pt','h_pt','s_pt','cv_pt','cp_pt','a_pt','jt_pt','mu_pt','r_pt','k_pt'}
        p=varargin{3}; 
        T=varargin{4}; 
        checkpT(p,T);
        p1=findp(p);
        
        switch prop(1)
            case 'v', val1=myi(dat.isoP{p1}.T,dat.isoP{p1}.v,T); val2=myi(dat.isoP{p1+1}.T,dat.isoP{p1+1}.v,T); %[m3/kg] Specific Volume
            case 'u', val1=myi(dat.isoP{p1}.T,dat.isoP{p1}.u,T); val2=myi(dat.isoP{p1+1}.T,dat.isoP{p1+1}.u,T); %[kJ/kg] Internal Energy
            case 'h', val1=myi(dat.isoP{p1}.T,dat.isoP{p1}.h,T); val2=myi(dat.isoP{p1+1}.T,dat.isoP{p1+1}.h,T); %[kJ/kg] Enthalpy
            case 's', val1=myi(dat.isoP{p1}.T,dat.isoP{p1}.s,T); val2=myi(dat.isoP{p1+1}.T,dat.isoP{p1+1}.s,T); %[kJ/kgK] Entropy
            case 'c'
                switch prop(2)
                    case 'v', val1=myi(dat.isoP{p1}.T,dat.isoP{p1}.cv,T); val2=myi(dat.isoP{p1+1}.T,dat.isoP{p1+1}.cv,T); %[kJ/kgK] Specific heat coefficient constant volume
                    case 'p', val1=myi(dat.isoP{p1}.T,dat.isoP{p1}.cp,T); val2=myi(dat.isoP{p1+1}.T,dat.isoP{p1+1}.cp,T); %[kJ/kgK] Specific heat coefficient constant pressure
                end
            case 'a', val1=myi(dat.isoP{p1}.T,dat.isoP{p1}.a,T); val2=myi(dat.isoP{p1+1}.T,dat.isoP{p1+1}.a,T); %[m/s] Sound speed
            case 'm', val1=myi(dat.isoP{p1}.T,dat.isoP{p1}.mu,T); val2=myi(dat.isoP{p1+1}.T,dat.isoP{p1+1}.mu,T); %[uPa*s] Viscosity
            case 'r', val1=myi(dat.isoP{p1}.T,dat.isoP{p1}.r,T); val2=myi(dat.isoP{p1+1}.T,dat.isoP{p1+1}.r,T); %[kg/m^3] Density
            case 'k', val1=myi(dat.isoP{p1}.T,dat.isoP{p1}.k,T); val2=myi(dat.isoP{p1+1}.T,dat.isoP{p1+1}.k,T); %[W/mK] Thermal Conductivity     
        end
        
        ret=interp1([dat.isoP{p1}.P(1) dat.isoP{p1+1}.P(1)],[val1 val2],p);
    case 't_ps'
        p=varargin{3}; 
        s=varargin{4};
        checkp(p);
        if (p<dat.Pcrit) 
            sl=INIST(varargin{1},'sl_p',p);
            sv=INIST(varargin{1},'sv_p',p);
            tsat=INIST(varargin{1},'tsat_p',p);
            if s>=sl && s<=sv 
                ret(1) = (s-sl)/(sv-sl);
                ret(2) = tsat;
            else
                eq=@(x) INIST(varargin{1},'s_pt',p,x)-s;
                options=optimset('Display','none');
                ret=fsolve(eq,tsat*1.1,options);
            end
        else
            eq=@(x) INIST(varargin{1},'s_pt',p,x)-s;
            options=optimset('Display','none');
            ret=fsolve(eq,dat.Tcrit*1.1,options);
        end
    case 't_ph'
        p=varargin{3}; 
        h=varargin{4};
        checkp(p);
        if (p<dat.Pcrit) 
            hl=INIST(varargin{1},'hl_p',p);
            hv=INIST(varargin{1},'hv_p',p);
            tsat=INIST(varargin{1},'tsat_p',p);
            if h>=hl && h<=hv 
                ret(1) = (h-hl)/(hv-hl);
                ret(2) = tsat;
            else
                eq=@(x) INIST(varargin{1},'h_pt',p,x)-h;
                options=optimset('Display','none');
                ret=fsolve(eq,tsat*1.1,options);
            end
        else
            eq=@(x) INIST(varargin{1},'h_pt',p,x)-h;
            options=optimset('Display','none');
            ret=fsolve(eq,dat.Tcrit*1.1,options);
        end
    case 'add_p'
        error('add_p is disabled')
        p=varargin{3};
        
        outr=0;
        try
            p1=findp(p);
        catch
            outr=1;
        end
        
        if outr==0 && ( ...
                abs(dat.isoP{p1}.P-p)<1e-3 || ...
                abs(dat.isoP{p1+1}.P-p)<1e-3 ) % don't add a pressure too close to a current isobar
            fprintf('p=%f is too close to a current isobar \n',p);
        else

            dat = get_isobar(dat, p, 8, 8); % download the isobar

            if length(varargin)>3
                    filename = varargin{4}; % Add the isobar to the database file if requested
                    file = load(filename); % Loads the requested database from the file
                    field = fieldnames(file); % Takes the name of the databse (Maybe it's not IND)

                    % At this point, dat contains the data specific to the
                    % compound, does not know anything about IND (E.g., IND.He)
                    % The isobaric data is added to this structure (E.g., IND.He)
                    % In ascending order of pressure
                    % Now the new data is added to the entire IND database.
                    % (field{:}) usually resolves to IND, but the name might
                    % be different in the future, hence this format.
                    % (dat.name) is the name of the compound, as it is in 
                    % the IND database (E.g., He)
                    file.(field{:}).(dat.name) = dat;

                    % And now the data loaded from file is saved onto a 
                    % Separate variable. This is needed because the save
                    % function of matlab does not like working with structure
                    % field names.
                    % This name, here IND, is the one that will be used when the
                    % structure is loaded once again in the future
                    IND = file.(field{:});

                    % Saves the database to the requested filename
                    % The content of the files should now be the previous
                    % database, as it was before, plus the new isobar
                    % NOTE: the name of the file can by anything, the
                    % load function will still load the structure
                    % as "IND". So, load('new_database.mat') will still
                    % load everything as "IND", because that's the second
                    % parameter of "save".
                    save(filename, 'IND');

                    % The following command loads the file to the base
                    % workspace, in order to be used in the next call
                    % evalin( 'base', 'load(filename)' );
            end
        end
        ret=dat;
    case 'minp'
        ret=dat.isoP{1}.P;
    case 'maxp'
        ret=dat.isoP{end}.P;
    case 'mint'
        ret=dat.isoP{1}.T(1);
    case 'maxt'
        ret=dat.isoP{1}.T(end);
    otherwise
        error('Unknown input parameter');
end


    function check_satp(p) 
        if p<dat.Psat(1)            
            error(sprintf('uhh? p=%e bar is too low, min sat. pressure for %s is %e bar',p,dat.name,dat.Psat(1)));
        end
        if (p>dat.Psat(end))
            error(sprintf('uhh? p=%e bar is above %s critical pressure=%e bar',p,dat.name,dat.Psat(end)));
        end
    end

    function check_undefined_p(p)
        if p==dat.Pcrit
            error(sprintf('You cannot get the values of cv, cp and a for critical pressure'));
        end
    end

    function check_undefined_t(t)
        if t==dat.Tcrit
            error(sprintf('You cannot get the values of cv, cp and a for critical temperature'));
        end
    end

    function check_satT(T)
        if T<dat.Tsat(1)
            error(sprintf('uhh? T=%e K is too low, min sat. temp. for %s is %e K',T,dat.name,dat.Tsat(1)));
        end
        if (T>dat.Tsat(end))
            error(sprintf('uhh? T=%e K is above %s critical temperature=%e K',T,dat.name,dat.Tsat(end)));
        end
    end

    function checkp(p)
        if ( (p < (dat.isoP{1}.P(1))) || (p > (dat.isoP{end}.P(1)) ) )
            error(sprintf('uhh? P=%e out of range of %s pressure (%e to %e bar)',p,dat.name,dat.isoP{1}.P,dat.isoP{end}.P));
%             error('uhh? P out of range');
        end
    end 

    function checkT(T)
        if ( (T < (dat.isoP{1}.T(1))) || (T > (dat.isoP{1}.T(end))) )
            error(sprintf('uhh? T=%e out of range of %s temperatures (%e to %e K)',T,dat.name,dat.isoP{1}.T(1),dat.isoP{1}.T(end)));
        end
    end 

    function checkpT(P,T) 
        checkp(P);
        checkT(T);
        % Now, we will check that P and T don't match saturation conditions
        if (P>dat.Psat(end)) % supercritical
            return;
        end
        if (T>dat.Tsat(end)) % supercritical
            return;
        end
        
        if ( T<=dat.Tsat(end) )
            Psat=INIST(varargin{1},'Psat_t',T);
            if Psat==P
                error('saturation !');
            end
        end
    end 

    function i=findp(p)
        for i=1:length(dat.isoP)-1
            if ( p >= dat.isoP{i}.P(1) ) && ( p <= dat.isoP{i+1}.P(1) )
                %fprintf('asked for p=%f, found between p=%f and p=%f \n',p,dat.isoP{i}.P,dat.isoP{i+1}.P);
                return
            end
        end
        error('Outside Splendiferous Pressure Range !');
    end

    function y=myi(X,Y,x) % interpolates avoiding l-v phase change singularity
        if x<X(1)
            error('Too low: Splendiferous error');
        end
        if x>X(end)
            error('Too high: Splendiferous error');
        end        
        q=find(diff(X)==0);
        if isempty(q) 
            y=interp1(X,Y,x);
            return
        end    
        if (x<=X(q))
            y=interp1(X(1:q),Y(1:q),x);
        else
            y=interp1(X(q+1:end),Y(q+1:end),x);
        end
    end

end

function [ fname ] = osi( fname )
% Manel Soria July 2019
% operating system independent
% given a path name
% changes / to \ or viceversa, only if needed, to suit the operating system
% eg osi('a/b/c') excuted in a windows machine will return 'a\b\c'

if ismac || isunix % to unix
    fname(fname=='\')='/';
else % windows
    fname(fname=='/')='\';    
end

end