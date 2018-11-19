%%% List of the model parameters %%%

Method = 'CV';          % FE - Finite Elements Method, CV - Control Volume Method for the solid

N = 100;                % number of nodes along X axis
M = 20;                 % number of nodes along radius of the particles
NL = 40;                % anode/separator boundary
NR = 60;                % separator/cathode boundary

L = 1e-4;               % cell thickness, m
L1 = 4e-5;              % anode thickness, m
Ldelta = 2e-5;          % separator thickness, m
L2 = L - L1 - Ldelta;   % cathode thickness, m
L1d = L1 + Ldelta;      % cathode boundary coordinate, m
a = 1e-5;               % electrode particle radius, m
A = 1e-4;               % electrode cross-sectional area, m^2

dx = 1/(N - 1);         % constant grid step along X axis (dimensionless)
dr = 1/(M - 1);         % constant grid step along radius of the particles (dimensionless)
hx = dx*L;              % constant grid step along X axis
hr = dr*a;              % constant grid step along radius of the particles

I0 = 1e-4;              % maximum current draw, A
t0 = 0;                 % charge/discharge period (set 0 if I = I0 = const), s
Rc = 3.58e-3;           % contact resistance, Ohm
tplus = 0.38;           % transference number

% Anode parameters
bet_a = pi/(2*a);       % BET (Brunauer-Emmett-Teller) surface area, m^-1
B_a = 0.3789;           % permeability factor of electrolyte
k0_a = 2e-11;           % reaction rate constant, m^5/2 mol^-1/2 s^-1
el_a = 0.5236;          % volume fraction of electrolyte
cmax_a = 20950;         % maximum concentration of Li ions in solid, mol m^-3
sigma_s_a = 5e-3;       % solid conductivity, S m^-1
Ds_a = 9e-16;           % solid ionic diffusivity, m^2 s^-1

% Cathode parameters
bet_c = pi/(2*a);       % BET (Brunauer-Emmett-Teller) surface area, m^-1
B_c = 0.3789;           % permeability factor of electrolyte
k0_c = 2e-11;           % reaction rate constant, m^5/2 mol^-1/2 s^-1
el_c = 0.5236;          % volume fraction of electrolyte
cmax_c = 20950;         % maximum concentration of Li ions in solid, mol m^-3
sigma_s_c = 5e-3;       % solid conductivity, S m^-1
Ds_c = 9e-16;           % solid ionic diffusivity, m^2 s^-1

% Separator parameters
B_s = 0.3789;           % permeability factor of electrolyte
el_s = 0.5236;          % volume fraction of electrolyte

% Solver parameters
V_stop = 2.0;           % stop simulation when the potential in solid (Cathode) reach Vstop, V
pre_Jacobian = false;   % use pre-generated Jacobian matrix (true/false)
rel_tol = 1e-5;         % relative tolerance of solver
abs_tol = 1e-5;         % absolute tolerance of solver
max_timestep = 10.0;    % maximum time step, s
t_end = 2.7007969e4;              % integrate equations from 0 to t_end (set 0 to perform a precursor calculation)

% Output parameters
n_lines = 14;           % number of lines to plot

% Initial conditions
c0_a = 1000;            % initial concentration of Li ions in the electrolyte in Anode, mol m^-3
cs0_a = 0.95*cmax_a;    % initial concentration of Li ions in solid in Anode, mol m^-3
c0_c = 1000;            % initial concentration of Li ions in the electrolyte in Cathode, mol m^-3
cs0_c = 0.05*cmax_c;    % initial concentration of Li ions in solid in Cathode, mol m^-3
c0_s = 1000;            % initial concentration of Li ions in the electrolyte in the separator, mol m^-3

% Other constants
F = 96487;              % Faraday's constant, C mol^-1
R = 8.3144;             % universal gas constant, J mol^-1 K^-1
T = 298;                % absolute temperature, K

% Derived indexes
Ndelta = NR - NL;                           % number of nodes for the separator
N_c_a = 3*N - Ndelta + 2;                   % starting index for concentration in solid in anode
N_c_c = 3*N - Ndelta + 2 + M*NL;            % starting index for concentration in solid in cathode
N_end = 2*N + (N - Ndelta + 1)*(M + 1);     % the last index

% Anode-Separator-Cathode combined parameters
bet = ASC_combine(bet_a, bet_c, 0, N, NL, NR);
B = ASC_combine(B_a, B_c, B_s, N, NL, NR);
k0 = ASC_combine(k0_a, k0_c, 0, N, NL, NR);
el = ASC_combine(el_a, el_c, el_s, N, NL, NR);
cmax = ASC_combine(cmax_a, cmax_c, 1, N, NL, NR);


function val = ASC_combine(val_a, val_c, val_s, N, NL, NR)
% Helper function to set Anode-Separator-Cathode combined parameters
    val = zeros(N, 1);
    val(1:NL-1) = val_a;
    val(NL:NR-1) = val_s;
    val(NR:N) = val_c;
end
