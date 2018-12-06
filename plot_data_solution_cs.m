%%% Plot concentration in solid %%%

disp('Plotting concentration in solid...')

figure('Name', 'C_s');

subplot(2, 3, 1);
plot(r.*a_a, y0(N_c_a:N_c_a+M-1), '.r');
hold on;
plot(r.*a_a, y(2:(n_lines-1), N_c_a:N_c_a+M-1), 'k');
hold on;
plot(r.*a_a, y(n_lines, N_c_a:N_c_a+M-1), '.b');
xlim([0 a_a]);
title('');
xlabel('r [m]');
ylabel('c_s(x=0) [mol/m^3]');
grid on;

subplot(2, 3, 2);
plot(r.*a_a, y0(N_c_a+M*fix(NL/2):N_c_a+M*fix(NL/2)+M-1), '.r');
hold on;
plot(r.*a_a, y(2:(n_lines-1), N_c_a+M*fix(NL/2):N_c_a+M*fix(NL/2)+M-1), 'k');
hold on;
plot(r.*a_a, y(n_lines, N_c_a+M*fix(NL/2):N_c_a+M*fix(NL/2)+M-1), '.b');
xlim([0 a_a]);
title('');
xlabel('r [m]');
ylabel('c_s(x=L_1/2) [mol/m^3]');
grid on;

subplot(2, 3, 3);
plot(r.*a_a, y0(N_c_a+M*(NL-1):N_c_a+M*NL-1), '.r');
hold on;
plot(r.*a_a, y(2:(n_lines-1), N_c_a+M*(NL-1):N_c_a+M*NL-1), 'k');
hold on;
plot(r.*a_a, y(n_lines, N_c_a+M*(NL-1):N_c_a+M*NL-1), '.b');
xlim([0 a_a]);
title('');
xlabel('r [m]');
ylabel('c_s(x=L_1) [mol/m^3]');
grid on;

subplot(2, 3, 4);
plot(r.*a_c, y0(N_c_c:N_c_c+M-1), '.r');
hold on;
plot(r.*a_c, y(2:(n_lines-1), N_c_c:N_c_c+M-1), 'k');
hold on;
plot(r.*a_c, y(n_lines, N_c_c:N_c_c+M-1), '.b');
xlim([0 a_c]);
title('');
xlabel('r [m]');
ylabel('c_s(x=L_1+\delta) [mol/m^3]');
grid on;

subplot(2, 3, 5);
plot(r.*a_c, y0(N_c_c+M*fix((N-NR)/2):N_c_c+M*fix((N-NR)/2)+M-1), '.r');
hold on;
plot(r.*a_c, y(2:(n_lines-1), N_c_c+M*fix((N-NR)/2):N_c_c+M*fix((N-NR)/2)+M-1), 'k');
hold on;
plot(r.*a_c, y(n_lines, N_c_c+M*fix((N-NR)/2):N_c_c+M*fix((N-NR)/2)+M-1), '.b');
xlim([0 a_c]);
title('');
xlabel('r [m]');
ylabel('c_s(x=L_1+\delta+L_2/2) [mol/m^3]');
grid on;

subplot(2, 3, 6);
plot(r.*a_c, y0(N_end-M+1:N_end), '.r');
hold on;
plot(r.*a_c, y(2:(n_lines-1), N_end-M+1:N_end), 'k');
hold on;
plot(r.*a_c, y(n_lines, N_end-M+1:N_end), '.b');
xlim([0 a_c]);
title('');
xlabel('r [m]');
ylabel('c_s(x=L) [mol/m^3]');
grid on;
