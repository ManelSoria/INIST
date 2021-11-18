% NFP - 
% Non-ideal Fluid Properties (previously INIST)
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% example: check NFP against XSteam, for water
% Requires XSteam, it can be downloaded from Mathworks File Exchange
% Remember to use addpath 

clear all
close all

tv=linspace(70,400,50); % temperatures (C) for each pressure
pv=[0.2,1,10,100,200]; % pressures

for j=1:length(pv)
        
        for i=1:length(tv)
            t=tv(i);
            sX(i)=XSteam('s_pt',pv(j),t);
            sI(i)=NFP('H2O','s_pt',pv(j),t+273.15);
        end

        plot(sX,tv,'r-','LineWidth',2);
        hold on
        plot(sI,tv,'bo');
    
end

ylabel('Temperature (oC)')
xlabel('Entropy (kJ/kgK)')
title('H2O Ts diagram for several pressures; XSteam (red) vs. NFP (blue) ')

set(gca,'FontSize',18)