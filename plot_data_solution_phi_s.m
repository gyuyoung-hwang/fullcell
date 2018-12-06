%%% Plot potential in solid %%%

disp('Plotting potential in solid...')

figure('Name', 'Potential');

clearvars y_axis
y_axis_a(1:NL, 1:n_lines) = 0;
for i = 1:n_lines
    for j = 1:NL
        y_axis_a(j, i) = y(i, 2*N+j);
    end
end

y_axis_c(1:N-NR+1, 1:n_lines) = 0;
for i = 1:n_lines
    for j = 1:N-NR+1
        y_axis_c(j, i) = y(i, 2*N+NL+j);
    end
end

subplot(1, 3, 1);
plot(t, y_axis_a(:, 1:n_lines), '-k');
xlim([0 t_max]);
title('');
xlabel('t [s]');
ylabel('\phi_s (Anode) [V]');
grid on;

subplot(1, 3, 2);
plot(t, y_axis_c(:, 1:n_lines), '-k');
xlim([0 t_max]);
title('');
xlabel('t [s]');
ylabel('\phi_s (Cathode) [V]');
grid on;

subplot(1, 3, 3);
plot(t, y(:, 3*N+NL-NR+1) - y(:, 2*N+1) + Rc*I_of_t(t, I0, t0), '-k');
xlim([0 t_max]);
title('');
xlabel('t [s]');
ylabel('Total voltage [V]');
grid on;

filename = 'fullcell_phi_s.xlsx';
warning('off', 'MATLAB:xlswrite:AddSheet');
xlswrite(filename, t, 'time');
xlswrite(filename, y_axis_a(:, 1:n_lines), 'phi_s (Anode)');
xlswrite(filename, y_axis_c(:, 1:n_lines), 'phi_s (Cathode)');
