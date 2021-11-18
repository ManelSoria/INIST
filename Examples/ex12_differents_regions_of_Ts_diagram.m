% NFP - 
% Non-ideal Fluid Properties (previously INIST)
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% Plot isobar 
close all


P1=1; T1=16;                    s1=NFP('H2','s_pt',P1,T1);
P2=1; T2=NFP('H2','tsat_p',P2); s2=NFP('H2','sl_p',P2);
P3=1; T3=NFP('H2','tsat_p',P2); s3=NFP('H2','sv_p',P3);
P4=1; T4=30;                    s4=NFP('H2','s_pt',P4,T4);

P5=40; T5=40;                   s5=NFP('H2','s_pt',P5,T5);
NFP_plotisobar('H2',P1,'b',2)
NFP_plotisobar('H2',P5,'k',2)

SZ=120;

hold on

scatter(s1,T1,SZ,'ob','filled')
scatter(s2,T2,SZ,'og','filled')
scatter(s3,T3,SZ,'ok','filled')
scatter(s4,T4,SZ,'om','filled')
scatter(s5,T5,SZ,'or','filled')
legend({'Subcritical isobar','Saturation bell','','Supercritical isobar','','','Subcooled liquid','Saturated liquid','Saturated vapor','Overheated vapor','Supercritical'},'Location','NW')
xlim([-10 30])
ylim([12 50])
set(gca,'FontSize',18)
