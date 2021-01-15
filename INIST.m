function [ret] = INIST(varargin)
% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% Units: T(K), p(bar), h and u: kJ/kg, v: m^3/kg, rho: kg/m^3 s: kJ/kgK,
% a: m/s, cv and cp: kJ/kgK, JT: bar/K, mu: Pa.s, k: W/mK, MM: kg/mol
% SF: N.m
% 1st argument: substance properties
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
% 
%   other functions:
%       't_hp', h, p, T0   returns the temperature given the enthalpy, pressure and a guessed temperature
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

try
    if isempty(IND) || ~isfield(IND,varargin{1})  
        set = load(varargin{1});
        IND.(varargin{1}) = set.(varargin{1});
    end
catch
    error('%s not found',varargin{1})
end


dat = IND.(varargin{1});
prop=lower(varargin{2}); 

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
        
%  computation of  temperature from enthalpy and pressure
    case 't_hp'
        species = varargin{1}; %species considered
        h = varargin{3}; %enthalpy [kJ/kg]
        p = varargin{4}; %pressure [bar]
        T_0 = varargin{5}; %temperature guess [K]
        
        options=optimset(...
        'Display','none',...
        'MaxIter',10000,...
        'TolFun', 1.0e-10,...
        'TolX',1.0e-4);

        P_crit=INIST(species,'pcrit'); %Critical pressure computation

        if p>P_crit
            has2b0 = @(T) h - INIST(species,'h_pt',p,T);
            ret = fsolve(has2b0,T_0,options);            
        else
            hl = INIST(species,'hl_p',p);
            hv = INIST(species,'hv_p',p);
            if h<=hv && h>=hl
                % In saturation line
                x2 = (h-hl)/(hv-hl);
                if x2>1 || x2<0 
                    error('Error in the vapour quality!!'); 
                else
                    ret = INIST(species,'tsat_p',p);
                end
            else
                % Not in the saturation line
                has2b0 = @(T) h - INIST(species,'h_pt',p,T);
                ret = fsolve(has2b0,T_0,options);
            end
        end
                
%  saturated liquid properties as a function of pressure    
    case 'vl_p',  p=varargin{3}; check_satp(p);ret=interp1(dat.Psat,dat.vl,p);        
    case 'ul_p',  p=varargin{3}; ccheck_satp(p);ret=interp1(dat.Psat,dat.ul,p);
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

        
    case {'v_pt','u_pt','h_pt','s_pt','cv_pt','cp_pt','a_pt','jt_pt','mu_pt','r_pt','k_pt'}
        p=varargin{3}; 
        T=varargin{4}; 
        checkp(p); 
        checkT(T); 
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
                ret=tsat;
            else
                eq=@(x) INIST(varargin{1},'s_pt',p,x)-s;
                options=optimset('Display','none');
                ret=fsolve(eq,tsat,options);
            end
        else
            eq=@(x) INIST(varargin{1},'s_pt',p,x)-s;
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
            error(sprintf('uhh? T=%e K is above %s critical pressure=%e K',T,dat.name,dat.Tsat(end)));
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

