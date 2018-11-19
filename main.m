%%% MAIN Program -- Entry point %%%

clear all
format long

% Physical and simulation parameters
parameters

if strcmp(Method, 'FE')
    disp('MAIN Program: Full Cell -- FE Method')
elseif strcmp(Method, 'CV')
    disp('MAIN Program: Full Cell -- FE+CV Method')
else
    disp('The Method should be either FE or CV in parameters.m')
end

% Non-dimensional computational grid
x = linspace(0, 1, N);
r = linspace(0, 1, M);

% Initial conditions
Us0_a = Ueq_of_cs_anode(cs0_a, cmax_a);         % initial potential in solid (Anode), V
Us0_c = Ueq_of_cs_cathode(cs0_c, cmax_c);       % initial potential in solid (Cathode), V
y0 = zeros(1, N_end);                           % reset all
y0(1:NL) = c0_a;                                % concentration of Li ions in liquid in Anode
y0(NL+1:NR-1) = c0_s;                           % concentration of Li ions in liquid in the separator
y0(NR:N) = c0_c;                                % concentration of Li ions in liquid in Cathode
y0(N+1:2*N) = 0;                                % potential in the electrolyte
y0(2*N+1:2*N+NL) = Us0_a;                       % potential in solid in Anode
y0(2*N+NL+1:3*N+1-Ndelta) = Us0_c;              % potential in solid in Cathode
y0(3*N+2-Ndelta:3*N+1-Ndelta+M*NL) = cs0_a;     % concentration of Li ions in solid in Anode
y0(3*N+2-Ndelta+M*NL:end) = cs0_c;              % concentration of Li ions in solid in Cathode

% Create options structure for the solver
disp('Initializing the solver...')
if t_end > 0
    if strcmp(Method, 'FE')
        opts = odeset('Mass',     @(t,y) mass_matrix_FEM(t, y, N, M, Ndelta, N_end), ...
                      'MaxStep',         max_timestep, ...
                      'RelTol',          rel_tol, ...
                      'AbsTol',          abs_tol, ...
                      'Stats',           'off');
    elseif strcmp(Method, 'CV')
        opts = odeset('Mass',     @(t,y) mass_matrix_FE_CV(t, y, N, M, Ndelta, N_end), ...
                      'MaxStep',         max_timestep, ...
                      'RelTol',          rel_tol, ...
                      'AbsTol',          abs_tol, ...
                      'Stats',           'off');
    end
else
    if strcmp(Method, 'FE')
        opts = odeset('Mass',     @(t,y) mass_matrix_FEM(t, y, N, M, Ndelta, N_end), ...
                      'Events',   @(t,y) auto_stop(t, y, N, Ndelta, V_stop, Rc, I0, t0), ...
                      'RelTol',          rel_tol, ...
                      'AbsTol',          abs_tol, ...
                      'Stats',           'off');
    elseif strcmp(Method, 'CV')
        opts = odeset('Mass',     @(t,y) mass_matrix_FE_CV(t, y, N, M, Ndelta, N_end), ...
                      'Events',   @(t,y) auto_stop(t, y, N, Ndelta, V_stop, Rc, I0, t0), ...
                      'RelTol',          rel_tol, ...
                      'AbsTol',          abs_tol, ...
                      'Stats',           'off');
    end
end

% Start precursor calculation if necessary
if t_end > 0
    t_max = t_end;
else
    disp('Starting precursor calculation...')
    tic
    t_out = linspace(0, 100000, 2);
    if strcmp(Method, 'FE')
        [t,y] = ode15s(@(t,y) scheme_FEM(t, y, x, r), t_out, y0, opts);
    elseif strcmp(Method, 'CV')
        [t,y] = ode15s(@(t,y) scheme_FE_CV(t, y, x, r), t_out, y0, opts);
    end
    time_pre = toc;
    t_max = max(t);
    disp(['Time: ' num2str(time_pre) ' seconds'])
    disp(['Discharge time: ' num2str(t_max) ' seconds (' num2str(t_max/3600) ' hours)'])
end

% Start calculation for the output
disp('Starting calculation for the output...')
tic
t_out = linspace(0, t_max, n_lines);
if strcmp(Method, 'FE')
    [t,y] = ode15s(@(t,y) scheme_FEM(t, y, x, r), t_out, y0, opts);
elseif strcmp(Method, 'CV')
    [t,y] = ode15s(@(t,y) scheme_FE_CV(t, y, x, r), t_out, y0, opts);
end
time_fin = toc;
disp(['Time: ' num2str(time_fin) ' seconds'])

% Output
n_lines = size(t,1);
plot_data_params
plot_data_solution
plot_data_solution_cs
plot_data_solution_phi_s
plot_data_control

disp('...All done')
disp('========================================')
