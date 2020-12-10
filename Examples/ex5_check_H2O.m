% INIST - 
% Interpolation of Nonideal Idiosyncratic Splendiferous Tables
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% example: check INIST against XSteam, for water
% Requires XSteam
clearvars
load('IND');
close
length(IND.H2O.isoP)

tv=linspace(70,400,50); % temperatures (C) for each pressure
pv=[0.2,1,10,100,200]; % pressures

for j=1:length(pv)
        %%%% NEED add_p
        %IND.H2O=INIST(IND.H2O,'add_p',pv(j)); % Download the pressure if not available
        %%%%
        
        for i=1:length(tv)
            t=tv(i);
            sX(i)=XSteam('s_pt',pv(j),t);
            sI(i)=INIST('H2O','s_pt',pv(j),t+273.15);
        end

        plot(sX,tv,'r-');
        hold on
        plot(sI,tv,'bo');
    
end

ylabel('Temperature (oC)')
xlabel('Entropy (kJ/kgK)')
title('H2O Ts diagram for several pressures; XSteam (red) vs. INIST (blue) ')

length(IND.H2O.isoP)

save('IND','IND');