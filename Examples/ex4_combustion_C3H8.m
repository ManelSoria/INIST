function ex4_combustion_C3H8
% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% Combustion LC3H8 / LO2
% Example: liquid phase propellants
% Requires HGS

species={'C3H8',...
    'CO2',...
    'CO',...
    'O2',...
    'O',...
    'H2',...
    'H',...
    'OH',...
    'H2O'};

% Propellants enters at the combustion chamber  as saturated liquid at
% chamber pressure
P0=20; % bar  

TO2=INIST('O2','tsat_P',P0)

TC3H8=INIST('C3H8','tsat_P',P0)

% 1 mol C3H8, 5 mol O2

np=[1;0;0;5;0;0;0;0;0];

% neq=hgseq(species,np,2000,P0)


deltah_O2=(INIST('O2','h_pT',P0,400)-INIST('O2','hl_P',P0))*INIST('O2','MM'); % kJ/mol
%deltah_O2=0 % uncomment for gas phase propellants
hO2 = HGSsingle('O2','h',400)-deltah_O2 % enthalpy at inlet, in hgs reference


deltah_C3H8=(INIST('C3H8','h_pT',P0,400)-INIST('C3H8','hl_P',P0))*INIST('C3H8','MM'); % kJ/mol
%deltah_C3H8=0  % uncomment for gas phase propellants
hC3H8 = HGSsingle('C3H8','h',400)-deltah_C3H8

Hin=np(1)*hO2+np(4)*hC3H8

[Tp,~,~,~] = HGStp(species,np,'H',Hin,P0)


end

