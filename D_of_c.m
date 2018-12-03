function D = D_of_c(c)
% Diffusion coefficient of electrolyte, m^2 s^-1

% Liquid ionic diffusivity (Valoen and Reimers)
%D = 5.253e-10*exp(-3.071e-4*c);  % m^2 s^-1

% Ecker, Käbitz, Laresgoiti et al.
D = 2.4e-10;
