clearvars

% From the paper
T=248.20; % K
p=703.10/1e5; % bar

rho=NFP('CO2','r_pt',p,T) % Density kg/m^3 
gamma=NFP('CO2','cp_pt',p,T)/NFP('CO2','cv_pt',p,T) % gamma 
mu=NFP('CO2','mu_pt',p,T) % dynamic viscosity 
Ru=0.0083144621; % kJ/(kgK)
Rg=1000*Ru/NFP('CO2','MM') % kJ/(KgK)