% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% plot C3H8 isobars

% check enthalpy difference in hgs and INIST

% both INIST and HGS have to be in the path 

clearvars
close



INIST_plotisobar('C3H8',[1,10])



% check that difference is similar, despite the reference value is
% different

HGSsingle('C3H8','h',400)- HGSsingle('C3H8','h',390) % kJ/mol
(INIST('C3H8','h_pt',1,400)-INIST('C3H8','h_pt',1,390))*INIST('C3H8','MM') % kJ/kg


% extend hgs to consider vaporisation enthalpy

% enthalpy difference between 400K at 1 bar and saturated liquid at 1 bar
deltaH=(INIST('C3H8','h_pt',1,400)-INIST('C3H8','hl_p',1))*INIST('C3H8','MM') % kJ/kg



