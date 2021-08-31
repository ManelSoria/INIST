% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% O2 turbopump work


clearvars


p1=4.48 % bar
p2=110.45 % bar
mdot=1807.7 % kg/s


T1=INIST('O2','Tsat_p',p1)
h1=INIST('O2','hl_p',p1)
s1=INIST('O2','sl_p',p1)


s2s=s1;
% T2s ?

T2s=INIST('O2','T_ps',p2,s2s)

h2s=INIST('O2','h_pt',p2,T2s)

wbs=h1-h2s % kJ/kgK

Wbs=mdot*wbs % kW




v1=INIST('O2','vl_p',p1) % m^3/kg

Wbs2=mdot*v1*(p1-p2)*1e5/1000 % kg/s * kg/m^3 * Pa = J /s -> kW
