% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% Creating IND database
% Please note that will require some hours and it do not save
% automatically. If your internet shut down, save IND manually and start
% from this point adjusting this code
% Regards, Caleb

species = {'H2O' 'N2' 'H2' 'CO' 'CO2' 'N2O' 'CH4O' 'CH4' 'C2H6' 'C2H4'...
    'C3H8' 'C3H6' 'C3H4' 'C4H10' 'C5H12' 'C6H14' 'C6H12' 'C6H6' 'NH3'...
    'He' 'O2' 'R134a'};
MM = [0.018 0.024 0.002 0.028 0.044 0.044 0.032 0.016 0.030 0.028 ...
    0.044 0.042 0.040 0.058 0.072 0.086 0.084 0.078 0.017 ...
    0.004 0.032 0.10203];
idcas = {'C7732185' 'C7727379' 'C1333740' 'C630080' 'C124389' 'C124389'...
    'C67561' 'C74828' 'C74840' 'C74851' 'C74986' 'C115071' 'C74997'...
    'C106978' 'C109660' 'C110543' 'C110827' 'C71432' 'C7664417'...
    'C7440597' 'C7782447' 'C811972'};


% This function takes as an input an existing database entry (e.g. IND.He)
% and  download the data directly from the NIST website

s_base_url = 'https://webbook.nist.gov/cgi/fluid.cgi?Action=Data';

digits = 8;

% Isobars
tinc = 2;
tmax = 5000;
tmin = 0;
s_Type = 'isoBar';

%New isobars
isobars1 = 0.001:0.001:0.01;
isobars2 = 0.02:0.01:0.1;
isobars3 = 0.2:0.1:10;
isobars4 = 10.5:0.5:100;
isobars5 = 101:1:250;
isobars6 = 255:5:500;
isobars = [isobars1 isobars2 isobars3 isobars4 isobars5 isobars6];

% Saturated
Tinc = 1;
Tmin = 0;
Tmax = 10000;
s_Type1 = 'SatP';

for ii=1:length(species)
    IND.(species{ii}).name = species{ii};
    IND.(species{ii}).idcas = idcas{ii};
    IND.(species{ii}).MM = MM(ii);
    fprintf('Downloading saturated data for %s... ',species{ii});
    % The webread function call the base url and adds the url-enconded 
    % parameters. The weboptions specifies a table stile for the return.
    satT = webread(s_base_url, ...
        'ID', idcas{ii}, 'Type', s_Type1, 'Digits', digits, 'TLow', Tmin, 'THigh', Tmax, 'TInc', Tinc, ...
        'RefState', 'DEF', 'TUnit', 'K', 'PUnit', 'bar', 'DUnit', 'kg/m3', 'HUnit', 'kJ/kg', ...
        'WUnit', 'm/s', 'VisUnit', 'Pa*s', 'STUnit', 'N/m', 'Wide', 'on', weboptions('ContentType','table'));
    
    % The table is then parsed
    IND.(species{ii}) = parseTableSaturated(IND.(species{ii}), satT);
    fprintf('OK\n');     
    
    for jj =1:length(isobars)
        fprintf('Downloading isobaric data for P = %.5f bar for %s... ', isobars(jj),species{ii});
        % The webread function call the base url and adds the url-enconded 
        % parameters. The weboptions specifies a table stile for the return
        isobar = webread(s_base_url, ...
            'ID', idcas{ii}, 'Type', s_Type, 'Digits', digits, 'P', isobars(jj), ...
            'RefState', 'DEF', 'TUnit', 'K', 'PUnit', 'bar', 'DUnit', 'kg/m3', 'HUnit', 'kJ/kg', ...
            'WUnit', 'm/s', 'VisUnit', 'Pa*s', 'STUnit', 'N/m', 'Wide', 'on', ...
            'TLow', tmin, 'THigh', tmax, 'TInc', tinc, weboptions('ContentType','table'));

        % The table is then parsed
        IND.(species{ii}) = parseTableIsobaric(IND.(species{ii}), jj, isobar);
        fprintf('OK\n');
    end
    
    % Eval is a dangerous function. Care!
    name = [species{ii} ' = IND.' species{ii} ';']; 
    eval(name);
    save(species{ii},species{ii});
end



%% Nested functions
function ret = parseTableIsobaric(ret, index, tab)

%     This function parses the table output of a webread for an isobar. 
%     It distributes the properties on the IND database using the 
%     appropriate field names, only if the column names
%     (hence the units) have not changed from the default
   
    for col = 1 : size(tab, 2)
        for row = 1 : size(tab, 1)
            switch tab.Properties.VariableNames{col}
                case 'Temperature_K_'
                    ret.isoP{index}.T(row) = new_s2d(tab.Temperature_K_(row));
                case 'Pressure_bar_'
                    ret.isoP{index}.P(1) = new_s2d(tab.Pressure_bar_(1));
                case 'Density_kg_m3_'
                    ret.isoP{index}.r(row) = new_s2d(tab.Density_kg_m3_(row));
                case 'Volume_m3_kg_'
                    ret.isoP{index}.v(row) = new_s2d(tab.Volume_m3_kg_(row));
                case 'InternalEnergy_kJ_kg_'
                    ret.isoP{index}.u(row) = new_s2d(tab.InternalEnergy_kJ_kg_(row));
                case 'Enthalpy_kJ_kg_'
                    ret.isoP{index}.h(row) = new_s2d(tab.Enthalpy_kJ_kg_(row));
                case 'Entropy_J_g_K_'
                    ret.isoP{index}.s(row) = new_s2d(tab.Entropy_J_g_K_(row));
                case 'Cv_J_g_K_'
                    ret.isoP{index}.cv(row) = new_s2d(tab.Cv_J_g_K_(row));
                case 'Cp_J_g_K_'
                    ret.isoP{index}.cp(row) = new_s2d(tab.Cp_J_g_K_(row));
                case 'SoundSpd__m_s_'
                    ret.isoP{index}.a(row) = new_s2d(tab.SoundSpd__m_s_(row));
                case 'Viscosity_Pa_s_'
                    ret.isoP{index}.mu(row) = new_s2d(tab.Viscosity_Pa_s_(row));
                case 'Therm_Cond__W_m_K_'
                    ret.isoP{index}.k(row) = new_s2d(tab.Therm_Cond__W_m_K_(row));
                case 'Phase'
                    ret.isoP{index}.phase(row) = new_s2d(tab.Phase(row));
            end
        end
    end  
end

function ret = parseTableSaturated(ret, tab)

%     This function parses the table output of a webread for the
%     saturated values. It distributes the properties on the IND
%     database using the appropriate field names, only if the column names
%     (hence the units) have not changed from the default
    
    for col = 1 : size(tab, 2)     
        for row = 1 : size(tab, 1)
            
            switch tab.Properties.VariableNames{col}
                case 'Temperature_K_'
                    ret.Tsat(row) = new_s2d(tab.Temperature_K_(row));
                    ret.Tcrit(1) = new_s2d(max(ret.Tsat));
                case 'Pressure_bar_'
                    ret.Psat(row) = new_s2d(tab.Pressure_bar_(row));
                    ret.Pcrit(1) = new_s2d(max(ret.Psat));
                case 'Density_l_Kg_m3_'
                    ret.rl(row) = new_s2d(tab.Density_l_Kg_m3_(row));
                case 'Density_v_Kg_m3_'
                    ret.rv(row) = new_s2d(tab.Density_v_Kg_m3_(row));
                case 'Volume_l_M3_kg_'
                    ret.vl(row) = new_s2d(tab.Volume_l_M3_kg_(row));
                case 'Volume_v_M3_kg_'
                    ret.vv(row) = new_s2d(tab.Volume_v_M3_kg_(row));
                case 'InternalEnergy_l_KJ_kg_'
                    ret.ul(row) = new_s2d(tab.InternalEnergy_l_KJ_kg_(row));
                case 'InternalEnergy_v_KJ_kg_'
                    ret.uv(row) = new_s2d(tab.InternalEnergy_v_KJ_kg_(row));
                case 'Enthalpy_l_KJ_kg_'
                    ret.hl(row) = new_s2d(tab.Enthalpy_l_KJ_kg_(row));
                case 'Enthalpy_v_KJ_kg_'
                    ret.hv(row) = new_s2d(tab.Enthalpy_v_KJ_kg_(row));
                case 'Entropy_l_J_g_K_'
                    ret.sl(row) = new_s2d(tab.Entropy_l_J_g_K_(row));
                case 'Entropy_v_J_g_K_'
                    ret.sv(row) = new_s2d(tab.Entropy_v_J_g_K_(row));
                case 'Cv_l_J_g_K_'
                    ret.cvl(row) = new_s2d(tab.Cv_l_J_g_K_(row));
                case 'Cv_v_J_g_K_'
                    ret.cvv(row) = new_s2d(tab.Cv_v_J_g_K_(row));
                case 'Cp_l_J_g_K_'
                    ret.cpl(row) = new_s2d(tab.Cp_l_J_g_K_(row));
                case 'Cp_v_J_g_K_'
                    ret.cpv(row) = new_s2d(tab.Cp_v_J_g_K_(row));
                case 'SoundSpd__l_M_s_'
                    ret.al(row) = new_s2d(tab.SoundSpd__l_M_s_(row));
                case 'SoundSpd__v_M_s_'
                    ret.av(row) = new_s2d(tab.SoundSpd__v_M_s_(row));
                case 'Viscosity_l_Pa_s_'
                    ret.mul(row) = new_s2d(tab.Viscosity_l_Pa_s_(row));
                case 'Viscosity_v_Pa_s_'
                    ret.muv(row) = new_s2d(tab.Viscosity_v_Pa_s_(row));
                case 'Therm_Cond__l_W_m_K_'
                    ret.kl(row) = new_s2d(tab.Therm_Cond__l_W_m_K_(row));
                case 'Therm_Cond__v_W_m_K_'
                    ret.kv(row) = new_s2d(tab.Therm_Cond__v_W_m_K_(row));
            end
        end 
    end
end


function output = new_s2d(input)
    if isstring(input)
        output = str2double(input);
    elseif iscell(input)
        output = str2double(input);
    else
        output = input;
    end
end
