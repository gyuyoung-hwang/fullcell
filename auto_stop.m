function [value, isterminal, direction] = auto_stop(t, y, N, Ndelta, V_stop, Rc, I0, t0)

% Simulation will stop when phi_s_cathode(x=L) == Vstop + Rc*I
value = y(3*N+1-Ndelta) - Rc*I_of_t(t, I0, t0) - V_stop;  % event function

% Simulation will stop when phi_s_anode(x=0) == Vstop + Rc*I
%value = y(2*N+1) - Rc*I_of_t(t, I0, t0) - V_stop;  % event function

isterminal = 1;  % terminate the integration if event function == 0
direction  = 0;

end
