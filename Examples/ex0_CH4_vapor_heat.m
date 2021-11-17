% NFP - 
% Non-ideal Fluid Properties (previously INIST)
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020


% Examples of NFP use with methane

clearvars 
% "clearvars" clears vars but not global vars (where NFP data is stored)
% you can use "clear all" instead but the data will be reloaded
% not a big deal 

close all

sp='CH4'; % species

% saturation temperature at 10 bar
NFP(sp,'Tsat_p',10) % K

% saturation pressure at 180 K
NFP(sp,'Psat_T',120) % bar

% saturation pressure at 300 K
NFP(sp,'Psat_T',180) % bar

% latent l->v heat a 20 bar
NFP(sp,'hv_p',20)-NFP(sp,'hl_p',20) % kJ/kg

% entropy of saturated liquid and vapour at 12 bar
NFP(sp,'sl_p',12) % kJ/kgK
NFP(sp,'sv_p',12) % kJ/kgK

% latent vaporisation heat vs. sat pressure
pv=linspace(1,NFP(sp,'pcrit'),50); % 50 values 
for i=1:length(pv)
 dhv(i)= NFP(sp,'hv_p',pv(i))-NFP(sp,'hl_p',pv(i));
end
figure

plot(pv,dhv,'LineWidth',2);
title('C3H8 vaporisation heat vs. saturation pressure');
xlabel('Psat (bar)');
ylabel('\Delta h (kJ/khK)');
grid
set(gca,'FontSize',18)

% Properties at 5 bar and 400 K
NFP(sp,'h_pt',5,400) % Enthalpy kJ/kgK
NFP(sp,'s_pt',5,400) % Entropy kJ/kgK
NFP(sp,'u_pt',5,400) % Internal energy kJ/kg
NFP(sp,'v_pt',5,400) % Volume m^3/kg

% CH4 is expanded isentropically from 100 bar, 200K to
% 5 bar. Determine final state and plot the expansion in a Ts diagram


NFP_plotisobar(sp,[ 1, 100],'k',2)
set(gca,'FontSize',18)


p1=100; % bar
T1=200; % K
p2=1; % bar
s1=NFP(sp,'s_pt',p1,T1) % kJ/kgK
s2=s1
sl=NFP(sp,'sl_p',p2)
sv=NFP(sp,'sv_p',p2)
x2=(s2-sl)/(sv-sl) % check that it is below 1 !
if x2>1
    error('huuu?? Not saturation conditions');
end
T2=NFP(sp,'Tsat_p',p2)
hold on
plot([s1 s2],[T1 T2],'-b','LineWidth',2);
SZ=120;
scatter(s1,T1,SZ,'og','filled')
scatter(s2,T2,SZ,'or','filled')

xlim([-1,5])
ylim([100,250])


 
