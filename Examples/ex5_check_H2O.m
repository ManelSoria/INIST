% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% example: check INIST against XSteam, for water
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
            sI(i)=INIST('H2O','s_pt',pv(j),t+273.15);
        end

        plot(sX,tv,'r-','LineWidth',2);
        hold on
        plot(sI,tv,'bo');
    
end

ylabel('Temperature (oC)')
xlabel('Entropy (kJ/kgK)')
title('H2O Ts diagram for several pressures; XSteam (red) vs. INIST (blue) ')

set(gca,'FontSize',18)