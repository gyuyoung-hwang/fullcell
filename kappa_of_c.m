function kappa = kappa_of_c(c)
% Ionic conductivity of electrolyte, S m^-1

% Liquid conductivity (Valoen and Reimers)
%kappa = 1e-4*c.*(5.2069096 - 0.002143628*c + 2.34402e-7*c.^2).^2;  % S m^-1

% Ecker, Käbitz, Laresgoiti et al.
cm = 1e-3*c;  % mol/m^3 to mol/l
kappa = 0.2667*cm.^3 - 1.2983*cm.^2 + 1.7919*cm + 0.1726;
