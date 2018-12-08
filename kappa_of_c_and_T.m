function kappa = kappa_of_c_and_T(c,T)
% Ionic conductivity of electrolyte, S m^-1

R = 8.3144; % Ideal gas constant

% Liquid conductivity (Valoen and Reimers)
%kappa = 1e-4*c.*(5.2069096 - 0.002143628*c + 2.34402e-7*c.^2).^2;  % S m^-1

% Ecker, Käbitz, Laresgoiti et al.
cm = 1e-3*c;  % mol/m^3 to mol/l
% First, calculate kappa(c,296)
kappa_c_296 = 0.2667*cm.^3 - 1.2983*cm.^2 + 1.7919*cm + 0.1726;
% Then add Arrhenius temperature dependence
C = 296*exp(17100/(R*296)); % Normalization constant
kappa = C*kappa_c_296*exp(-17100/(R*T))./T;