function ex10_H2_storage
% HYDROGEN STORAGE
% Manel Soria 2021

% The quantity of H2 that a car would need to have a range of 300 miles
% has been estimated to be between 5 and 7 kg
% (https://doi.org/10.1016/j.ijhydene.2019.07.044)

% A tank contains 7 kg of H2 at -252.87°C and 1 atm (state 1)

% 1-State of the fluid

clear 
close all

p1=1.01325;  % bar
T1=-252.87+273.15; % K
m1=7; % kg
Tsat1=INIST('H2','tsat_p',p1);
if T1<Tsat1 
    fprintf('State1 is subcooled liquid (below the saturation temperature for its pressure\n');
else
    error('uhh?? the rest of the code assumes state 1 is Subcooled');
end

% 2-H2 density and volume of the tank

rho1=INIST('H2','r_pT',p1,T1); % kg/m^3
V1=m1/rho1; % m^3

fprintf('State1 density is %.1f kg/m^3 and tank volume %.1f l\n',rho1,V1*1000);

% NOTE: IN PRACTICE, SOME GASEOUS H2 WOULD BE PRESENT AT THE TOP OF THE
% TANK, THE PURPOSE OF THIS EXAMPLE IS TO SHOW THAT A CRYOGENIC LIQUID IS
% DIFFICULT TO STORE

% 3-Compare the tank volume with the volume that would be needed to store
% the same mass of liquid water (assume 1000 kg/m^3 density of water)

fprintf('It has to be %.1f times bigger than an equivalent water tank \n',... 
    V1*1000 /m1);

% 4-Assume that the insulation of the tank is such that 100W/m^2 of heat
% from the ambient enter it. Calculate the pressure of the tank as a
% function of time. Assume the tank is a perfectly rigid sphere.

% (The 100W/m^2 estimation for a LH2 tank is from 
% Thomas D. Bostock Ralph G. Scurlock
% Low-Loss Storage and Handling of Cryogenic Liquids 
% The Application of Cryogenic Fluid Dynamics
% Another source  doi:10.1088/1757-899X/171/1/012063) gives a higher value
% 200W/m^2
% Of course it depends on the insulation thickness

% V= (4/3)*pi*r^3
r=(V1 / ( (4/3)*pi )) ^(1/3); % tank radius (m)

S=4*pi*r^2; % tank surface m^2

q=S*100; % W heat flux

% If the tank is rigid, the 1st Thermodynamic principle reads:
% U2-U1 = Q 

U1=m1*INIST('H2','u_pt',p1,T1); % kJ

dt=10;

Q=dt*q; % 

U2=U1+Q; 
V2=V1;
m2=m1;

u2=U2/m2;
v2=V2/m2;

% So, given u2 and v2, we have to find p2 and T2:

options1=optimset('tolfun',1e-11,'display','none'); % a lower tolerance has to be set as v values are small

% FIRST APPROACH:

p2=fsolve(@htbz2,30,options1);  % find p2

htbz1=@(T2) INIST('H2','v_pt',p2,T2)-v2; % Given T2 and the current value of p2, htbz1 returns the error in specific volume

T2=fsolve(htbz1,30,options1); % recover T2

    function err=htbz2(p2) % given a pressure p2, finds the temperature T2 so that v=v2 and then returns the error in u

        % find T2 that has v2 at the given p2
        % The function has to be defined here because of the use of p2
        htbz1=@(T2) INIST('H2','v_pt',p2,T2)-v2;
        
        T2=fsolve(htbz1,30,options1);
        
        % find the difference to the internal energy 
        err=INIST('H2','u_pt',p2,T2)-u2; 
    end

% SECOND APPROACH: SOLVE BOTH EQUATIONS TOGETHER
% This is the most natural method but be carefully with the factor needed
% in the specific volume equation

    function errB=htbzB(x)
        p=x(1);
        T=x(2);
        % Both equations errors have to be roughly of the same order of magnitude
        % So (after some tests) we multiply the first by 1e4
        errB(1)=1e4*(INIST('H2','v_pt',p,T)-v2); 
        errB(2)=INIST('H2','u_pt',p,T)-u2;
    end

    xs=fsolve(@htbzB,[40,40],options1);

    p2_prime=xs(1);
    T2_prime=xs(2);

    fprintf('p2=%e bar T2=%e K \n',p2_prime,T2_prime);

    fprintf('Error in v = %e \n',INIST('H2','v_pt',p2,T2)-v2 );
    fprintf('Error in u = %e \n',INIST('H2','u_pt',p2,T2)-u2 );

    fprintf("difference in p2=%e \n",p2-p2_prime);
    fprintf("difference in T2=%e \n",T2-T2_prime);


end