% ORTHO-PARA HYDROGEN

clear 
close all

% The book "Ignition!: An Informal History of Liquid Rocket Propellants" 
% by  J.D.Clark says:
% "Each mole of hydrogen (2 grams) which changed from the ortho to the 
% para state gave off 337 calories of heat in the process. And since it 
% takes only 219 calories to vaporize one mole of hydrogen, you were in 
% real trouble. For if you liquefied a mass of hydrogen, getting a liquid 
% that was still almost three quarters ortho- hydrogen, the heat of the 
% subsequent transition of that to para-hydro-gen was enough to change the
% whole lot right back to the gaseous state. All without the help of any 
% heat leaking in from the outside

% Verify these claims using INIST and the comment in the book
% Thermodynamic properties of cryogenic fluids by Jacobsen, 
% Richard T. Leachman, Jacob W. Lemmon, Eric W. Penoncello, Steven G 
% "For all comparisons and mixing with parahydrogen, the reference enthalpy 
% and entropy of saturated liquid orthohydrogen at the normal boiling point
% should be changed to 702.98 kJ/kg and 0.018269 kJ/kg-K, respectively."

% First we find the saturation temperatures at 1atm
P=1.013; % bar
T1p=INIST('pH2','tsat_p',P)
T1o=INIST('H2','tsat_p',P)

% Now, we check the enthalpy values for both substances as sat liquid, 1atm
H1p=INIST('pH2','hl_p',P) % kJ/kg
H1o=INIST('H2','hl_p',P) % kJ/kg

% We see that both are close to zero, while according to J.D.Clark the
% difference should be significant. The answer of course is that for each
% substance the reference enthalpy is different.

% So, we impose
H1o=702.98;

% And now, the transition isobaric heat is:
DeltaHop=H1p-H1o

% The heat is released in the transition from ortho (obtained from the gas
% state at 300K) to para, the stable form at 20K. 75% of the gas is ortho.

% So, for 1mol, (about 2g), the heat released (in absolute value) is
Qop=abs(2e-3*(H1p-H1o)*0.239*1000) % cal

% Very close to the J.D.Clark value

% Considering that only 3/4 of the gas is ortho, the transition heat is
Qop=0.75*Qop 

% Next check is the vaporization heat of 2 g for the para form:

Qvap=2e-3*(INIST('pH2','hv_p',P)-INIST('pH2','hl_p',P) )*0.239*1000

% also very close 

% Meaning that ALL the H2 just liquified would evaporate even if the
% liquid could be stored in a perfectly adiabatic container
