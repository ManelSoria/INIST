% NFP - 
% Non-ideal Fluid Properties (previously INIST)
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% Example: R134a vapour compression refrigeration
%
%
%* 3.- Un cicle de refrigeracio per compressio de vapor funciona amb 
% 0.1 kg/s de R134a. El punt 1, situat a la sortida de l?evaporador, 
% es trova en estat de vapor saturat a la pressi? de 0.24 MPa . 
% Es suposa que el compressor ?s isentr?pic. 
% A la sortida del compressor (punt 2), la pressi? ?s de 0.7 MPa. 
% Al condensador la pressi? es mant? pr?cticament constant, i a la sortida 
% (punt 3) l?estat ?s de l?quid saturat. La v?lvula d?expansi? ?s isoent?lpica. 
% A l?evaporador la pressi? tamb? es mant? pr?cticament constant. 
% Es demana:
% a-Calcular el treball de compressi?
% b-Calcular el t?tol a la sortida de la v?lvula d?expansi?.
% c-Calcular la calor que despren el condensador.
% d-Calcular la calor que absorbeix l?evaporador.
% e-Calcular el COP

error('This example has to be reviewed');

clearvars
close all

plow=2.4; % bar
phigh=7;
mp=0.1; % kg/s

NFP_plotisobar('R134a',[plow,phigh]);

p1=plow;
T1=NFP('R134a','tsat_p',p1);
s1=NFP('R134a','sv_p',p1);
h1=NFP('R134a','hv_p',p1);
plot(s1,T1,'or');

p2=phigh;
s2=s1;
T2=NFP('R134a','t_ps',p2,s2);
h2=NFP('R134a','h_pt',p2,T2);
plot(s2,T2,'or');

p3=phigh;
T3=NFP('R134a','tsat_p',p3);
s3=NFP('R134a','sl_p',p3);
h3=NFP('R134a','hl_p',p3);
plot(s3,T3,'or');

p4=plow;
h4=h3;
T4=NFP('R134a','tsat_p',p4);
hl=NFP('R134a','hl_p',p4);
hv=NFP('R134a','hv_p',p4);
sl=NFP('R134a','sl_p',p4);
sv=NFP('R134a','sv_p',p4);
x4=(h4-hl)/(hv-hl);
s4=sl+x4*(sv-sl);
plot(s4,T4,'or');

axis([1 1.8 260 310])

wc=mp*(h1-h2) % kW
qc=mp*(h3-h2) % kW
qf=mp*(h1-h4) % kW

cop=abs(qf/wc);

fprintf('wc=%f kW\n',wc);
fprintf('x4=%f \n',x4);
fprintf('qc=%f kW \n',qc);
fprintf('qf=%f kW \n',qf);
fprintf('suma(q)-w = %f\n',qc+qf-wc)
fprintf('cop=%f \n',cop);

