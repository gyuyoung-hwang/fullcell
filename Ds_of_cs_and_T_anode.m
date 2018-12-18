function Ds = Ds_of_cs_and_T_anode(cs, cmax, T)
% Diffusion coefficient in solid for anode, m^2/s

R = 8.3144; % universal gas constant, J mol^-1 K^-1

x = cs/cmax;

% Graphite anode from Ecker, Kabitz, Laresgoiti et al. (GITT data, not EIS)
% Analytical fit (WebPlotDigitizer + MATLAB curve fiting tool)

% First, calculate kappa(cs ,cmax, 296) using exponential fit to EIS data
Ds_cs_296 = 2.275e-13*exp(-8.588*x);
% Now add Arrhenius temperature dependence
Ds = Ds_cs_296*exp(-30300/(R*T))*exp(30300/(R*296));

% Constant for testing purpose
%Ds = 1e-14;