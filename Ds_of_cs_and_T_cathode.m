function Ds = Ds_of_cs_and_T_cathode(cs, cmax, T)
% Diffusion coefficient in solid for cathode, m^2/s

R = 8.3144; % universal gas constant, J mol^-1 K^-1

x = cs/cmax;

% NMC406 cathode from Ecker, Kabitz, Laresgoiti et al.
% Assume no cs dependence in the first instance

% First, calculate kappa(cs ,cmax, 296) using exponential fit to EIS data
Ds_cs_296 = 6.3e-13;
% Now add Arrhenius temperature dependence
Ds = Ds_cs_296*exp(-80600/(R*T))*exp(80600/(R*296));