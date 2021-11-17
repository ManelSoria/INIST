% NFP - 
% Non-ideal Fluid Properties (previously INIST)
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% Plot isobar 
close all
NFP_plotisobar('H2',1,'b')
NFP_plotisobar('H2',30,'k')
xlim([-10 30])
ylim([0 50])
hold on
sliquid = NFP('H2','s_pt',1,16)
svapor = NFP('H2','s_pt',1,40)
sliquidsat = NFP('H2','sl_p',1)
svaporsat = NFP('H2','sv_p',1)
tsat = NFP('H2','tsat_p',1)
scatter(sliquid,16,50,'oc','filled')
scatter(sliquidsat,tsat,50,'og','filled')
scatter(svaporsat,tsat,50,'ok','filled')
scatter(svapor,40,50,'om','filled')
legend('1 bar isobar','Saturation bell','','30 bar isobar','','','Liquid T=12','Liquid saturated','Vapor saturated','Vapor T=40')
