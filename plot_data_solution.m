%%% Plot solution data %%%

disp('Plotting solution data...')

xa = linspace(0, L1, NL);
xc = linspace(L1d, L, N - NR + 1);

figure('Name', 'Solution');

subplot(3, 2, 1);
plot(x.*L, y0(1:N), '.r');
hold on;
plot(x.*L, y(2:(n_lines-1), 1:N), 'k');
hold on;
plot(x.*L, y(n_lines, 1:N), '.b');
xlim([0 L]);
title('');
xlabel('x [m]');
ylabel('c [mol/m^3]');
grid on;

subplot(3, 2, 2);
plot(x.*L, y0(N+1:2*N), '.r');
hold on;
plot(x.*L, y(2:(n_lines-1), N+1:2*N), 'k');
hold on;
plot(x.*L, y(n_lines, N+1:2*N), '.b');
xlim([0 L]);
title('');
xlabel('x [m]');
ylabel('\phi [V]');
grid on;

clearvars y_axis_0 y_axis
for i = 1:NL
    y_axis_0(i) = y0(N_c_a + i*M - 1);
    for j = 1:n_lines
        y_axis(j, i) = y(j, N_c_a + i*M - 1);
    end
end
subplot(3, 2, 3);
plot(xa, y_axis_0, '.r');
hold on;
plot(xa, y_axis(2:(n_lines-1), :), 'k');
hold on;
plot(xa, y_axis(n_lines, :), '.b');
xlim([0 L1]);
xlabel('x [m]');
ylabel('c_s(r=a) (Anode) [mol/m^3]');
grid on;

clearvars y_axis_0 y_axis
for i = 1:N-NR+1
    y_axis_0(i) = y0(N_c_c + i*M - 1);
    for j = 1:n_lines
        y_axis(j, i) = y(j, N_c_c + i*M - 1);
    end
end
subplot(3, 2, 4);
plot(xc, y_axis_0, '.r');
hold on;
plot(xc, y_axis(2:(n_lines-1), :), 'k');
hold on;
plot(xc, y_axis(n_lines, :), '.b');
xlim([L1d L]);
xlabel('x [m]');
ylabel('c_s(r=a) (Cathode) [mol/m^3]');
grid on;

subplot(3, 2, 5);
plot(xa, y0(2*N+1:2*N+NL), '.r');
hold on;
plot(xa, y(2:(n_lines-1), 2*N+1:2*N+NL), 'k');
hold on;
plot(xa, y(n_lines, 2*N+1:2*N+NL), '.b');
xlim([0 L1]);
xlabel('x [m]');
ylabel('\phi_s (Anode) [V]');
grid on;

subplot(3, 2, 6);
plot(xc, y0(2*N+NL+1:3*N+1-Ndelta), '.r');
hold on;
plot(xc, y(2:(n_lines-1), 2*N+NL+1:3*N+1-Ndelta), 'k');
hold on;
plot(xc, y(n_lines, 2*N+NL+1:3*N+1-Ndelta), '.b');
xlim([L1d L]);
xlabel('x [m]');
ylabel('\phi_s (Cathode) [V]');
grid on;
