function D = D_of_c_and_T(c,T)
F = 96487; % Faraday's constant, C/mol
kB = 1.381e-23; % Boltzmann's constant, J/K
q_e = 1.602e-19; % Elementary charge, C
% Liquid ionic diffusivity (Ecker et al. 2015)
D = (kB/(F*q_e))*kappa_of_c_and_T(c,T).*T./c;  % m^2 s^-1
