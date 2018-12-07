%%% Plot V/Q characteristic %%%

disp('Plotting V/Q characteristic...')

figure('Name', 'V/Q characteristic');

V_anode(1:n_lines) = 0;
for i = 1:n_lines
    V_anode(i) = y(i, 2*N+1);
end

V_cathode(1:n_lines) = 0;
for i = 1:n_lines
    V_cathode(i) = y(i, 3*N+NL-NR+1);
end

plot(48*t.*I_of_t(t, I0, t0)/3600, V_cathode-V_anode+Rc*I_of_t(t, I0, t0), '-k');
title('');
xlabel('Capacity [A h]');
ylabel('Total voltage [V]');
grid on;