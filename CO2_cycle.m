%Example of CO2 cooling cycle by Jordi Luque Barcons 2021
clc;
clear all;
close all;

%Given the cooling cycle shown in the figure problem_data.PNG, find the
%unknowns (if any) at each station of the cycle and plot the enthalpy
%variation in a h-P diagram:

%Problem to solve:
problem=imread('CO2_problem_data.PNG');
imshow(problem);

%Temperature and pressure data
T1=116.2+273.15; %Compressor exit temperature [K]
P1=105.06; %Compressor exit pressure [bar]
T2=100.71+273.15; %Condenser inlet temperature [K]
P2=105.06; %Condenser inlet pressure [bar]
T3=48.22+273.15; %Condenser exit temperature [K]
P3=104.4; %Condenser exit pressure [bar]
T4=47.69+273.15; %Expansion valve inlet temperature [K]
P4=104.4; %Expansion valve inlet pressure [bar]
P5=26.44; %Expansion valve exit pressure [bar]
T6=-9.88+273.15; %Evaporator inlet temperature [K]
P6=26.44; %Evaporator inlet pressure [bar]
T7=-5.37+273.15; %Evaporator exit temperature [K]
P7=25.8; %Evaporator exit pressure [bar]
T8=-2.39+273.15; %Compressor inlet temperature [K]
P8=25.8; %Compressor inlet pressure [bar]

%Draw the h-p diagram of the cycle:

% Properties at the exit of the compressor:
x = saturated_check(T1,P1,'CO2');

if x == -1
    h1=INIST('CO2','h_pt',P1,T1);
else
    hl1=INIST('CO2', 'hl_p',P1);
    hv1=INIST('CO2','hv_p',P1);
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the entrance of the condenser:
x = saturated_check(T2,P2,'CO2');

if x == -1
    h2=INIST('CO2','h_pt',P2,T2);
else
    hl2=INIST('CO2', 'hl_p',P2);
    hv2=INIST('CO2','hv_p',P2);
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the exit of the condenser:
x = saturated_check(T3,P3,'CO2');

if x == -1
    h3=INIST('CO2','h_pt',P3,T3);
else
    hl3=INIST('CO2', 'hl_p',P3);
    hv3=INIST('CO2','hv_p',P3);
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the entrance of the isoenthalpic valve:
x = saturated_check(T4,P4,'CO2');

if x == -1
    h4=INIST('CO2','h_pt',P4,T4);
else
    hl4=INIST('CO2', 'hl_p',P4);
    hv4=INIST('CO2','hv_p',P4);
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the exit of the isoenthalpic valve:
h5=h4; %Enthalpy is conserved
%We check if inside the saturation bell
hl5=INIST('CO2', 'hl_p',P5);
hv5=INIST('CO2', 'hv_p',P5);

if hl5<h5 && h5>hv5
    x=1; %Inside saturation bell
    T=INIST('CO2','tcrit');
else
    x=-1; %Outside saturation bell
    options=optimset(...
        'Display','none',...
        'MaxIter',10000,...
        'TolFun', 1.0e-10,...
        'TolX',1.0e-4);
    has2b0 = @(T) h5 - INIST('CO2','h_pt',P5,T);
    T5 = fsolve(has2b0,T4,options);
end
    
% Properties at the inlet of the evaporator
x = saturated_check(T6,P6,'CO2');

if x == -1
    h6=INIST('CO2','h_pt',P6,T6);
else
    hl6=INIST('CO2', 'hl_p',P6);
    hv6=INIST('CO2','hv_p',P6);
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the exit of the evaporator
x = saturated_check(T7,P7,'CO2');

if x == -1
    h7=INIST('CO2','h_pt',P7,T7);
else
    hl7=INIST('CO2', 'hl_p',P7);
    hv7=INIST('CO2','hv_p',P7);
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the inlet of the compressor
x = saturated_check(T8,P8,'CO2');

if x == -1
    h8=INIST('CO2','h_pt',P8,T8);
else
    hl8=INIST('CO2', 'hl_p',P8);
    hv8=INIST('CO2','hv_p',P8);
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

%Plot the enthalpy variation
P=[P1,P2,P3,P4,P5,P6,P7,P8,P1];
H=[h1,h2,h3,h4,h5,h6,h7,h8, h1];
figure
plot(H,P,'-o');
hold on;
labels = cellstr(num2str([1:length(P)-1]'));
text(H(1:end-1),P(1:end-1),labels,'VerticalAlignment','bottom','HorizontalAlignment','right');

%Plot saturation bell
Pcrit=INIST('CO2','pcrit');
p=linspace(6,Pcrit,1000);
for i=1:length(p)
    hl(i)=INIST('CO2', 'hl_p',p(i));
    hv(i)=INIST('CO2', 'hv_p',p(i));
end
plot(hl,p,'k');
hold on
plot(hv,p,'k');
grid minor
xlabel('h [kJ/kg]','Fontsize',14);
ylabel('P [bar]','Fontsize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nested function
 
function [x] = saturated_check(T,P,species)
% Function to check if a substance is inside the saturation bell or not
% Input data::
%   T --> Temperature [K]
%   P --> Pressure [bar]
%   rho --> substance density [kg/m^3]
%   species --> type of element [string]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output data::
%   x --> quality factor (between 0 and 1, or -1 if is not saturated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_crit = INIST(species,'tcrit');
P_crit = INIST(species,'pcrit');

%If the pressure and/or temperature are above the critical one, it is not
%saturated
if (P > P_crit || T > T_crit) 
    x = -1; 
    return;
end

%If the temperature is not equal to the saturated
%temperature for that given pressure it is not saturated
T_sat = INIST(species,'tsat_p',P);
if (T~=T_sat)
    x = -1;
    return;
end

x = 1; %Otherwise, inside saturation bell

end