% NFP - 
% Non-ideal Fluid Properties (previously INIST)
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% Plot isobar 


clearvars
close all

isobarV=[4 6 8 10 ];
NFP_plotisobar('C3H8',isobarV);
 set(gca,'FontSize',18)

 % add plot legend
LL={};
for i=1:numel(isobarV)
    LL{end+1}=sprintf('%.1f bar',isobarV(i));
end

legend(LL)
