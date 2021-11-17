% NFP - 
% Non-ideal Fluid Properties (previously INIST)
% (c) Manel Soria, Caleb Fuster, Lorenzo Frezza
% Data downloaded from NIST web page
% ESEIAAT - UPC - 2014-2020
%
% Example by  Alejandro Nunez Labielle and Mario Gayete


% A valve is connected to a O2 tank. 
% In the oxygen tank, we have saturated liquid at a pressure equal to
% 10 bar. At the valve outlet, the pressure has decreased to 7 bar. 
% We are askeed to obtain the O2 state after the valve, assuming that it is
% isoenthalpic

clearvars
close all
      
        
%Calcular les propietats al estat inicial 
p1=10;%bar pressio
h1=NFP('O2','hl_p',p1);%kJ/kg
s1=NFP('O2','sl_p',p1);%kJ/kgK

%Calcul estat 2
p2=7;%pressio en bar
h2=h1;
hl2=NFP('O2','hl_p',p2);
hv2=NFP('O2','hv_p',p2);
x2=(h2-hl2)/(hv2-hl2);
if x2>1 || x2<0
    Error('ugg. This title is not valid ');
end
% Check that entropy increases
T2=NFP('O2','tsat_p',p2);
sl2=NFP('O2','sl_p',p2);
sv2=NFP('O2','sv_p',p2);
s2=sl2+x2*(sv2-sl2);
Inc_s=s2-s1;
if Inc_s<0
    Error('ugg. Increment of entropy negative');
end
fprintf('Saturated vapour at T2=%f K \n',T2);
fprintf('Outlet quality is %8.4f\n',x2);
fprintf('Delta s = %8.4f kJ/kgK \n',s2-s1);
    

% Diagram
pv = linspace(1,NFP('O2','pcrit'),40); %bar

vl_p = NFP('O2','vl_p',pv);
vv_p = NFP('O2','vv_p',pv);

v1 = NFP('O2','vl_p',p1);


% v2 value (inside the saturation bell !)
v2 = NFP('O2','vl_p',p2) + x2*(NFP('O2','vv_p',p2)-NFP('O2','vl_p',p2));

fprintf('Specific volume at the exit = %8.4f m^3/kg \n',v2);

figure;
semilogx(vl_p,pv,'r','Linewidth',2);
hold on;
text(v1*1.05,p1+2,'1','Fontsize',14);
text(v2*1.05,p2+2,'2','Fontsize',14);
semilogx(vv_p,pv,'r','Linewidth',2);
semilogx(v1,p1,'ob','Linewidth',2);
semilogx(v2,p2,'ob','Linewidth',2);
semilogx([v1 v2],[p1 p2],'k','Linewidth',2);
xlabel('v','Fontsize',14);
ylabel('p','Fontsize',14);
set(gca,'Fontsize',12)