% NFP - 
% Non-ideal Fluid Properties (previously INIST)
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2021
%
% Creating NFP database
% Please note that will require some hours to download all and it do not 
% save automatically. If your internet shut down, save IND manually and
% start from this point adjusting this code
% Regards, Caleb


function Download_NFPdata(varargin)
% MM for all and reduced case
H_MM = 1.00794 / 1000;
O_MM = 15.9994 / 1000;
N_MM = 14.0067 / 1000;
C_MM = 12.0107 / 1000;
He_MM = 4.002602 / 1000;
F_MM = 18.9984 / 1000;

% Ref for all and reduced case
Ref1 = 'NBP state ref. Sets the enthalpy and entropy to zero for the saturated liquid at the normal boiling point temperature.';
Ref2 = 'IIR state ref. Sets to 200 kJ/kg and 1 kJ/(kg-K) for enthalpy and entropy, respectively, for the saturated liquid at 0°C';
Ref3 = 'Reference sets to H = 2551.013479 kJ/kg and S = 9.103679 J/g*K  at 300.0 K and 0.010 bar.';
Ref4 = 'Reference sets to H = 309.269908 kJ/kg and S = 6.839209 J/g*K  at 298.15 K and 1.000 bar.';
Ref5 = 'Reference sets to H = 1699.663687 kJ/kg and S = 9.356091 J/g*K  at 300.0 K and 0.010 bar.';
Ref6 = 'Reference sets to H = 271.013 kJ/kg andS = 6.41058 J/g*K  at 298.15 K and 1.000 bar.';

% isobars (bar)
isobars1 = 0.001:0.001:0.01;
isobars2 = 0.02:0.01:0.1;
isobars3 = 0.2:0.1:10;
isobars4 = 11:1:100;
isobars5 = 105:5:250;
isobars6 = 260:10:290;
isobars7 = 300:100:700;

iso1 = [isobars1 isobars2 isobars3 isobars4 isobars5 isobars6 isobars7];
iso2 = [isobars1 isobars2 isobars3 isobars4 isobars5];
    
if isempty(varargin) || (ischar(varargin{1}) && strcmpi(varargin{1},'all'))
    
    % Stopper in case of you missed something
    sure = input('ARE YOU SURE TO GO ON DOWNLOADING ALL DATA ?? \n');
    if ~strcmpi(sure,'Yes')
        return
    end
    
    species =  {'H2O' 'N2' 'H2' 'CO' 'CO2' 'N2O' 'CH4O' 'CH4' 'C2H6'...
                'C2H4' 'C3H8' 'C3H6' 'C3H4' 'C4H10' 'C5H12' 'C6H14'...
                'C6H12' 'C6H6' 'NH3' 'He' 'O2' 'R134a' 'pH2'};
    
    % Shout-out to Alessandro Rampazzo
    %              H2O        %N2       %H2        CO       CO2       N2O
    idcas = {'C7732185' 'C7727379' 'C1333740' 'C630080' 'C124389' 'C10024972'...
        ...%CH4O     CH4    C2H6     C2H4     C3H8      C3H6     C3H4
        'C67561' 'C74828' 'C74840' 'C74851' 'C74986' 'C115071' 'C74997'...
        ...%C4H10   C5H12     C6H14     C6H12     C6H6       NH3
        'C106978' 'C109660' 'C110543' 'C110827' 'C71432'  'C7664417'...
        ...%He         O2     R134a(C2H2F4)  pH2 (para-H2)
        'C7440597' 'C7782447'   'C811972'     'B5000001'}; 
    

    
    %calculating every MM using element MMs
    %        H2O         N2      H2         CO       CO2           N2O
    MM = [H_MM*2+O_MM, N_MM*2, H_MM*2, C_MM+O_MM, C_MM+O_MM*2, N_MM*2+O_MM, ...
        ...%  CH4O              CH4          C2H6           C2H4
        C_MM+H_MM*4 + O_MM,  C_MM+H_MM*4, C_MM*2+H_MM*6, C_MM*2+H_MM*4,...
        ...%  C3H8          C3H6         C3H4           C4H10
        C_MM*3+H_MM*8, C_MM*3+H_MM*6, C_MM*3+H_MM*4, C_MM*4+H_MM*10,...
        ...%  C5H12            C6H14         C6H12             C6H6
        C_MM*5+H_MM*12, C_MM*6 + H_MM*14, C_MM*6+H_MM*12, C_MM*6+H_MM*6,...
        ...%  NH3     He     O2          R134a(C2H2F4)   pH2 (para-H2)
        N_MM+H_MM*3, He_MM, O_MM*2, C_MM*2+H_MM*2+F_MM*4  H_MM*2];
    
          %H2O   %N2   %H2   CO    CO2    N2O  CH4O   CH4   C2H6   C2H4
    Ref ={Ref3  Ref4  Ref1  Ref1  Ref2   Ref1  Ref1   Ref1  Ref1   Ref1 ...
        ...%  C3H8   C3H6    C3H4   C4H10  C5H12  C6H14  C6H12  C6H6
              Ref2   Ref2    Ref2    Ref2   Ref1   Ref1  Ref1   Ref1  ...
        ...%  NH3     He     O2   R134a(C2H2F4)  pH2 (para-H2)
              Ref5   Ref1   Ref6      Ref2        Ref1 };
          
               %H2O   %N2   %H2    CO    CO2   N2O  CH4O   CH4   C2H6  C2H4
    isobars = { iso1, iso1, iso1, iso1, iso1, iso2, iso1, iso1, iso1, iso1, ...
      ...%  C3H8   C3H6    C3H4   C4H10  C5H12  C6H14  C6H12  C6H6
           iso1,   iso1,   iso2,  iso2,  iso2,   iso2,  iso2, iso2,  ...
      ...%  NH3     He     O2   R134a(C2H2F4)  pH2 (para-H2) 
           iso1,   iso1,  iso1,     iso1,        iso1 }; 
     
elseif (ischar(varargin{1}) && strcmpi(varargin{1},'reduced'))
    % Stopper in case of you missed something
    sure = input('ARE YOU SURE TO GO ON DOWNLOADING reduced DATA ?? \n');
    if ~strcmpi(sure,'Yes')
        return
    end
    
    species =  { 'H2O' 'N2' 'H2' 'CO2' 'CH4' 'He' 'O2' 'R134a'};
    
    % Shout-out to Alessandro Rampazzo
    %              H2O        %N2       %H2       CO2
    idcas = { 'C7732185' 'C7727379' 'C1333740' 'C124389' ...
         ... CH4      He          O2      R134a(C2H2F4)   
           'C74828' 'C7440597' 'C7782447'  'C811972'}; ...
    
    %calculating every MM using element MMs
    %       H2O         N2      H2          CO2        
    MM = [H_MM*2+O_MM, N_MM*2, H_MM*2,  C_MM+O_MM*2,  ...
        ...%    CH4     He       O2      R134a(C2H2F4)
          C_MM+H_MM*4, He_MM, O_MM*2, C_MM*2+H_MM*2+F_MM*4];  
    
          %H2O   %N2   %H2    CO2    CH4  
    Ref ={Ref3  Ref4  Ref1   Ref2    Ref1   ...
        ...%  He     O2   R134a(C2H2F4)
              Ref1   Ref6      Ref2};
          
                %H2O   %N2   %H2   CO2    CH4  
    isobars = { iso1, iso1, iso1,  iso1,  iso1,  ...
              ...%   He     O2   R134a(C2H2F4)     
                    iso1,  iso1,    iso1}; % 
elseif (ischar(varargin{1}) && strcmpi(varargin{1},'H2-O2'))
    % Stopper in case of you missed something
    sure = input('ARE YOU SURE TO GO ON DOWNLOADING H2 and O2 DATA ?? \n');
    if ~strcmpi(sure,'Yes')
        return
    end
    
    species =  {'H2' 'O2'};
    
    % Shout-out to Alessandro Rampazzo
    %            %H2        O2 
    idcas = { 'C1333740' 'C7782447'}; ...
    
    %calculating every MM using element MMs
    %      H2      % O2      
    MM = [H_MM*2, O_MM*2];  
    
          %H2    O2
    Ref ={Ref1  Ref6};
          
                %H2    O2 
    isobars = { iso1, iso1}; 
elseif (ischar(varargin{1}) && strcmpi(varargin{1},'H2-pH2'))
    % Stopper in case of you missed something
    sure = input('ARE YOU SURE TO GO ON DOWNLOADING H2 and pH2 DATA ?? \n');
    if ~strcmpi(sure,'Yes')
        return
    end
    
    species =  {'H2' 'pH2'};
    
    % Shout-out to Alessandro Rampazzo
    %            %H2        pH2 
    idcas = { 'C1333740' 'B5000001'}; ...
    
    %calculating every MM using element MMs
    %      H2      % pH2      
    MM = [H_MM*2, H_MM*2];  
    
          %H2    pH2
    Ref ={Ref1  Ref1};
          
                %H2    pH2 
    isobars = { iso1, iso1};    
elseif (ischar(varargin{1}) && strcmpi(varargin{1},'Own_data'))
    species = varargin{2};
    MM = varargin{3};
    idcas = varargin{4};
    Ref = varargin{5};
    isobars = varargin{6};
else
    error('You passed a data download case that it is not supported')
end

% Base link
s_base_url = 'https://webbook.nist.gov/cgi/fluid.cgi?Action=Data';

% Number digits
digits = 8;

% Isobaric parse
tinc = 0;
tmax = 5000;
tmin = 0;
s_Type = 'isoBar';

% Saturated parse
Pinc = 0;
Plow = 0;
PHigh = 10000;
s_Type1 = 'SatT';

for ii=1:length(species)
    
    fprintf('Downloading basic and saturated data for %s... ',species{ii});    
    
    IND.(species{ii}).name = species{ii};
    IND.(species{ii}).idcas = idcas{ii};
    IND.(species{ii}).MM = MM(ii);
    IND.(species{ii}).Ref = Ref{ii};
    fprintf('OK\n');   
    
    fprintf('Downloading saturated data... ');    
    % The webread function call the base url and adds the url-enconded 
    % parameters. The weboptions specifies a table stile for the return.
    satP = webread(s_base_url, ...
        'ID', idcas{ii}, 'Type', s_Type1, 'Digits', digits, 'Plow', Plow, 'PHigh', PHigh, 'PInc', Pinc, ...
        'RefState', 'DEF', 'TUnit', 'K', 'PUnit', 'bar', 'DUnit', 'kg/m3', 'HUnit', 'kJ/kg', ...
        'WUnit', 'm/s', 'VisUnit', 'Pa*s', 'STUnit', 'N/m', 'Wide', 'on', weboptions('ContentType','table','Timeout',15));

    % The table is then parsed
    IND.(species{ii}) = parseTableSaturated(IND.(species{ii}), satP);
    fprintf('OK\n');  
    
    for jj =1:length(isobars{ii})
        fprintf('Downloading isobaric data for P = %.3f bar for %s... ', isobars{ii}(jj),species{ii});
        % The webread function call the base url and adds the url-enconded 
        % parameters. The weboptions specifies a table stile for the return
        isobar = webread(s_base_url, ...
            'ID', idcas{ii}, 'Type', s_Type, 'Digits', digits, 'P', isobars{ii}(jj), ...
            'RefState', 'DEF', 'TUnit', 'K', 'PUnit', 'bar', 'DUnit', 'kg/m3', 'HUnit', 'kJ/kg', ...
            'WUnit', 'm/s', 'VisUnit', 'Pa*s', 'STUnit', 'N/m', 'Wide', 'on', ...
            'TLow', tmin, 'THigh', tmax, 'TInc', tinc, weboptions('ContentType','table','Timeout',15));

        % The table is then parsed
        IND.(species{ii}) = parseTableIsobaric(IND.(species{ii}), jj, isobar);
        fprintf('OK\n');
    end
    
   % Not a very good practice. Care!
    name = [species{ii} ' = IND.' species{ii} ';']; 
    eval(name);
    save(species{ii},species{ii});
end
end

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
                    ret.isoP{index}.phase{row} = new_s2d(tab.Phase(row));
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
        output = input{1};
        if ~isnan(str2double(output))
            output = str2double(output);
        elseif strcmp(output,'undefined')
            output = NaN;
        end
    else
        output = input;
    end
end