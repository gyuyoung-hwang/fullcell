%%% Plot model functions %%%

pts = 1000;  % number of points to plot (100 by default)

disp('Plotting model functions...')

figure('Name', 'Model Functions');

x_axis = linspace(0, 1, pts);
y_axis = Ueq_of_cs_anode(x_axis.*cmax_a, cmax_a);
subplot(2, 3, 1);
plot(x_axis, y_axis, '.k');
title('');
xlabel('c_s/c_{max}');
ylabel('U_{eq} (Anode) [V]');
grid on;

x_axis = linspace(0, 2*c0_a, pts);
y_axis = D_of_c(x_axis);
subplot(2, 3, 2);
plot(x_axis, y_axis, '.k');
ylim([1e-10 5e-10]);
title('');
xlabel('c [mol/m^3]');
ylabel('D (Anode) [m^2/s]');
grid on;

x_axis = linspace(0, 2*c0_a, pts);
y_axis = kappa_of_c(x_axis);
subplot(2, 3, 3);
plot(x_axis, y_axis, '.k');
title('');
xlabel('c [mol/m^3]');
ylabel('\kappa (Anode) [S/m]');
grid on;

x_axis = linspace(0, 1, pts);
y_axis = Ueq_of_cs_cathode(x_axis.*cmax_c, cmax_c);
subplot(2, 3, 4);
plot(x_axis, y_axis, '.k');
%ylim([0 8]);
title('');
xlabel('c_s/c_{max}');
ylabel('U_{eq} (Cathode) [V]');
grid on;

x_axis = linspace(0, 2*c0_c, pts);
y_axis = D_of_c(x_axis);
subplot(2, 3, 5);
plot(x_axis, y_axis, '.k');
ylim([1e-10 5e-10]);
title('');
xlabel('c [mol/m^3]');
ylabel('D (Cathode) [m^2/s]');
grid on;

x_axis = linspace(0, 2*c0_c, pts);
y_axis = kappa_of_c(x_axis);
subplot(2, 3, 6);
plot(x_axis, y_axis, '.k');
title('');
xlabel('c [mol/m^3]');
ylabel('\kappa (Cathode) [S/m]');
grid on;
