%Example of CO2 cooling cycle by Jordi Luque Barcons January 2021
clc;
clear all;
close all;

%Given the CO2 cooling cycle shown in the figure problem_data.PNG, find 
%the enthalpy at each point of the cooling cycle and plot it in a h-P
%diagram along with the saturation bell. The valve can be modelled as an isoenthalpic one.
%Additionally, find the temperature at the exit of the expansion valve. For
%a more detailed explanation of the cooling cycle see:

%https://www.overleaf.com/2284234274ptchfdjhyvqk

%Problem to solve:
problem=imread('CO2_problem_data.PNG'); %Extracted from:
%Joaquim Rigola, "Aplicació de Matlab-Octave a Problemes d'Enginyeria
%Tèrmica - Class notes, Escola Tècnica Superior d'Enginyeria Industrial,
%Aeroespacial i Audiovisual de Terrassa, UPC, Abril 2019

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
T6=-8.88+273.15; %Evaporator inlet temperature [K]
P6=26.44; %Evaporator inlet pressure [bar]
T7=-5.37+273.15; %Evaporator exit temperature [K]
P7=25.8; %Evaporator exit pressure [bar]
T8=-2.39+273.15; %Compressor inlet temperature [K]
P8=25.8; %Compressor inlet pressure [bar]

%Enthalpy at each point of the cycle

% Properties at the exit of the compressor:
x = diff_sat_T(T1,P1,'CO2');
if x == -1
    h1=INIST('CO2','h_pt',P1,T1);
else
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the entrance of the condenser:
x = diff_sat_T(T2,P2,'CO2');
if x == -1
    h2=INIST('CO2','h_pt',P2,T2);
else
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the exit of the condenser:
x = diff_sat_T(T3,P3,'CO2');
if x == -1
    h3=INIST('CO2','h_pt',P3,T3);
else
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the entrance of the isoenthalpic valve:
x = diff_sat_T(T4,P4,'CO2');
if x == -1
    h4=INIST('CO2','h_pt',P4,T4);
else
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the exit of the isoenthalpic valve:
h5=h4; %Enthalpy is conserved
%We check if inside the saturation bell
hl5=INIST('CO2', 'hl_p',P5); %Enthalpy of saturated liquid
hv5=INIST('CO2', 'hv_p',P5); %Enthalpy of saturated vapour

if hl5<h5 && h5<hv5
    x=1; %Inside saturation bell
    T5=INIST('CO2','tsat_p',P5);
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
x = diff_sat_T(T6,P6,'CO2');
if x == -1
    h6=INIST('CO2','h_pt',P6,T6);
else
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the exit of the evaporator
x = diff_sat_T(T7,P7,'CO2');
if x == -1
    h7=INIST('CO2','h_pt',P7,T7);
else
    error('Watch out! You need the vapour quality to find the enthalpy!');
end

% Properties at the inlet of the compressor
x = diff_sat_T(T8,P8,'CO2');
if x == -1
    h8=INIST('CO2','h_pt',P8,T8);
else
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
 
function [x] = diff_sat_T(T,P,species)
% Function to check if a substance is inside the saturation bell or not
% Input data::
%   T --> Temperature [K]
%   P --> Pressure [bar]
%   species --> type of element [string]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output data::
%   x --> boolean -1 for unsaturated 1 for saturated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_crit = INIST(species,'tcrit');
P_crit = INIST(species,'pcrit');

%If the pressure and/or temperature are above the critical one, it is not
%saturated
if (P > P_crit || T > T_crit) 
    x = -1; 
    return;
end

T_sat = INIST(species,'tsat_p',P); %Saturation temperature at that specific pressure

%If the difference between the input temperature and the saturation
%temperature is below 0.5K, saturated conditions will be considered

if abs(T-T_sat)<0.5 
    x = 1;
    return;
end

x = -1; %Otherwise, outside saturation bell

end