%%% Plot solution control data %%%

disp('Plotting solution control data...')

figure('Name', 'Control');

clearvars y_axis
y_axis(1:n_lines) = 0;
for i = 1:n_lines
    for j = 1:N
        y_axis(i) = y_axis(i) + y(i, j)*dx*el(j);
    end
end
subplot(1, 2, 1);
plot(t, 100*(y_axis-y_axis(1))/y_axis(1), '-or', 'MarkerFaceColor', 'r');
xlim([0 t_max]);
title('Lithium Conservation in liquid');
xlabel('Time [s]');
ylabel('Percentage deviation from initial value [%]');
grid on;

y_axis(1:n_lines) = 0;
for i = 1:n_lines
    for k = 1:NL
        indx = N_c_a+(k-1)*M;
        y_axis(i) = y_axis(i) + y(i, indx)*dr^3*dx*pi/6;
        for j = 2:M-1
            indx = N_c_a+(k-1)*M-1+j;
            y_axis(i) = y_axis(i) + y(i, indx)*pi/6*((r(j) + r(j+1))^3 - (r(j-1) + r(j))^3)*dx;
        end
        indx = N_c_a+(k-1)*M-1+M;
        y_axis(i) = y_axis(i) + y(i, indx)*4*pi/3*(1 - (1 - dr/2)^3)*dx;
    end
    for k = NR:N
        indx = N_c_c+(k-NR)*M;
        y_axis(i) = y_axis(i) + y(i, indx)*dr^3*dx*pi/6;
        for j = 2:M-1
            indx = N_c_c+(k-NR)*M-1+j;
            y_axis(i) = y_axis(i) + y(i, indx)*pi/6*((r(j) + r(j+1))^3 - (r(j-1) + r(j))^3)*dx;
        end
        indx = N_c_c+(k-NR)*M-1+M;
        y_axis(i) = y_axis(i) + y(i, indx)*4*pi/3*(1 - (1 - dr/2)^3)*dx;
    end
end
subplot(1, 2, 2);
plot(t, 100*(y_axis-y_axis(1))/y_axis(1), '-or', 'MarkerFaceColor', 'r');
xlim([0 t_max]);
title('Lithium Conservation in solid');
xlabel('Time [s]');
ylabel('Percentage deviation from initial value [%]');
grid on;
