function h=ex14_elevation(Tc)
% NFP - 
% Non-ideal Fluid Properties (previously INIST)
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020

% Obtain the elevation from the boiling point of water in ºC
p=NFP('H2O','psat_t',Tc+273.15);

options=optimset('display','none');
h=fsolve(@hastobezero,1000,options);

function r=hastobezero(h)
    [~,~,P,~]=atmosisa(h);
    P=P/1e5;
    r=P-p;
end
end
