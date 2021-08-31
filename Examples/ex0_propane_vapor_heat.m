% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020


% Examples of INIST use with propane

clearvars 
% "clearvars" clears vars but not global vars (where INIST data is stored)
% you can use "clear all" instead but the data will be reloaded
% not a big deal 

close all

 % saturation temperature at 10 bar
 INIST('C3H8','Tsat_p',10) % K
 
 % saturation pressure at 200 K
 INIST('C3H8','Psat_T',200) % bar
 
 % saturation pressure at 300 K
 INIST('C3H8','Psat_T',300) % bar

 % latent l->v heat a 20 bar
 INIST('C3H8','hv_p',20)-INIST('C3H8','hl_p',20) % kJ/kg
 
 % entropy of saturated liquid and vapour at 12 bar
 INIST('C3H8','sl_p',12) % kJ/kgK
 INIST('C3H8','sv_p',12) % kJ/kgK
 
 % latent vaporisation heat vs. sat pressure
 pv=linspace(1,INIST('C3H8','pcrit'),40); % 40 values 
 for i=1:length(pv)
     dhv(i)= INIST('C3H8','hv_p',pv(i))-INIST('C3H8','hl_p',pv(i));
 end
 plot(pv,dhv,'LineWidth',2);
 title('C3H8 vaporisation heat vs. saturation pressure');
 xlabel('Psat (bar)');
 ylabel('\Delta h (kJ/khK)');
 grid
 set(gca,'FontSize',18)
 
 % Properties at 5 bar and 400 K
INIST('C3H8','h_pt',5,400) % Enthalpy kJ/kgK
INIST('C3H8','s_pt',5,400) % Entropy kJ/kgK
INIST('C3H8','u_pt',5,400) % Internal energy kJ/kg
INIST('C3H8','v_pt',5,400) % Volume m^3/kg

% C3H8 is expanded isentropically from 5 bar, 280K to
% 1 bar. Determine final state
p1=5; % bar
T1=280; % K
p2=1; % bar
s1=INIST('C3H8','s_pt',p1,T1) % kJ/kgK
s2=s1
sl=INIST('C3H8','sl_p',p2)
sv=INIST('C3H8','sv_p',p2)
x2=(s2-sl)/(sv-sl) % check that it is below 1 !
if x2>1
    error('huuu?? Not saturation conditions');
end
T2=INIST('C3H8','Tsat_p',p2)



 
