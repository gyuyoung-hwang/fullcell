function kappa = kappa_of_c(c)
% Liquid conductivity (Valoen and Reimers)
kappa = 1e-4*c.*(5.2069096 - 0.002143628*c + 2.34402e-7*c.^2).^2;  % S m^-1
