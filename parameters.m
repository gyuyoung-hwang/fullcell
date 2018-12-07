%%% List of the model parameters %%%

Method = 'FE';          % FE - Finite Elements Method, CV - Control Volume Method for the solid

Nmm = 1;                % number of nodes per 10^-6 m in liquid
M = 20;                 % number of nodes along radius of the particles

L = 148.0e-6;           % cell thickness, m
L1 = 74.0e-6;           % anode thickness, m
Ldelta = 19.0e-6;       % separator thickness, m
L2 = L - L1 - Ldelta;   % cathode thickness, m
L1d = L1 + Ldelta;      % cathode boundary coordinate, m

N = round(L*1e6*Nmm) + 1;               % number of nodes along X axis
NL = round(L1*1e6*Nmm) + 1;             % anode/separator boundary
NR = round((L1 + Ldelta)*1e6*Nmm) + 1;  % separator/cathode boundary

dx = 1/(N - 1);         % constant grid step along X axis (dimensionless)
dr = 1/(M - 1);         % constant grid step along radius of the particles (dimensionless)
hx = dx*L;              % constant grid step along X axis, m

I0 = 0.15625;           % maximum current draw (1C current = 0.15625 A), A
t0 = 0;                 % charge/discharge period (set 0 if I = I0 = const), s
Rc = 0;                 % total contact resistance, Ohm
tplus = 0.26;           % transference number
C_SEI = 0.068;          % capacity loss due to SEI

% Anode parameters
a_a = 13.7e-6;          % anode particle radius, m
act_a = 1 - 0.445;      % active part of anode
A_a = 8.585e-3;         % anode cross-sectional area, m^2
hr_a = dr*a_a;          % constant grid step along radius of the particles, m
el_a = 0.329;           % volume fraction of electrolyte
bet_a = 3*(1-el_a)/a_a; % BET (Brunauer-Emmett-Teller) surface area, m^-1
k0_a = 1.45e-11;        % reaction rate constant, m^5/2 mol^-1/2 s^-1
B_a = el_a/2.03;        % permeability factor of electrolyte
cmax_a = 31920*act_a;   % maximum concentration of Li ions in solid, mol m^-3
sigma_s_a = 14;         % solid conductivity, S m^-1
Ds_a = 1.0e-14;         % solid ionic diffusivity, m^2 s^-1 -- should depend on c

% Cathode parameters
a_c = 6.49e-6;          % cathode particle radius, m
act_c = 1 - 0.42;       % active part of cathode
A_c = 8.585e-3;         % cathode cross-sectional area, m^2
hr_c = dr*a_c;          % constant grid step along radius of the particles, m
el_c = 0.296;           % volume fraction of electrolyte
bet_c = 3*(1-el_c)/a_c; % BET (Brunauer-Emmett-Teller) surface area, m^-1
k0_c = 3.01e-11;        % reaction rate constant, m^5/2 mol^-1/2 s^-1
B_c = el_c/1.94;        % permeability factor of electrolyte
cmax_c = 48580*act_c;   % maximum concentration of Li ions in solid, mol m^-3
sigma_s_c = 68.1;       % solid conductivity, S m^-1
Ds_c = 1.0e-13;         % solid ionic diffusivity, m^2 s^-1 -- should depend on c

% Separator parameters
A_s = A_a;              % separator cross-sectional area, m^2
el_s = 0.508;           % volume fraction of electrolyte
B_s = el_s/1.67;        % permeability factor of electrolyte

% Solver parameters
V_stop = 1.0;           % stop simulation when the potential in solid (Anode) reach Vstop, V
pre_Jacobian = false;   % use pre-generated Jacobian matrix (true/false)
rel_tol = 1e-5;         % relative tolerance of solver
abs_tol = 1e-5;         % absolute tolerance of solver
max_timestep = 10.0;    % maximum time step, s
t_end = 0;              % integrate equations from 0 to t_end (set 0 to perform a precursor calculation)

% Output parameters
n_lines = 21;           % number of lines to plot

% Initial conditions
c0_a = 1000;            % initial concentration of Li ions in the electrolyte in Anode, mol m^-3
                        % initial concentration of Li ions in solid in Anode, mol m^-3
cs0_a = 0.60*cmax_c*L2*(1-el_c)*act_c/(L1*(1-el_a)*act_a);
%cs0_a = (0.74-C_SEI)*cmax_c*L2*(1-el_c)/(L1*(1-el_a));
c0_c = 1000;            % initial concentration of Li ions in the electrolyte in Cathode, mol m^-3
cs0_c = 0.26*cmax_c;    % initial concentration of Li ions in solid in Cathode, mol m^-3
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
