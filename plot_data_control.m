%%% Plot solution control data %%%

disp('Plotting solution control data...')

figure('Name', 'Control');

y_axis_liq(1:n_lines) = 0;
for i = 1:n_lines
    for j = 1:NL
        y_axis_liq(i) = y_axis_liq(i) + y(i, j)*hx*el_a*A_a;
    end
    for j = NL+1:NR-1
        y_axis_liq(i) = y_axis_liq(i) + y(i, j)*hx*el_s*A_s;
    end
    for j = NR:N
        y_axis_liq(i) = y_axis_liq(i) + y(i, j)*hx*el_c*A_c;
    end
end
subplot(1, 3, 1);
plot(t, 100*(y_axis_liq - y_axis_liq(1))/y_axis_liq(1), '-or', 'MarkerFaceColor', 'b');
%plot(t, y_axis_liq, '-or', 'MarkerFaceColor', 'b');  % absolute values
xlim([0 t_max]);
title('Lithium Conservation in liquid');
xlabel('Time [s]');
ylabel('Percentage deviation from initial value [%]');
grid on;

y_axis_s(1:n_lines) = 0;
for i = 1:n_lines
    for k = 1:NL
        indx = N_c_a+(k-1)*M;
        y_axis_s(i) = y_axis_s(i) + y(i, indx)*dr^3/6*bet_a*a_a*A_a/4*hx;
        for j = 2:M-1
            indx = N_c_a+(k-1)*M-1+j;
            y_axis_s(i) = y_axis_s(i) + y(i, indx)/6*((r(j) + r(j+1))^3 - (r(j-1) + r(j))^3)*bet_a*a_a*A_a/4*hx;
        end
        indx = N_c_a+(k-1)*M-1+M;
        y_axis_s(i) = y_axis_s(i) + y(i, indx)*4/3*(1 - (1 - dr/2)^3)*bet_a*a_a*A_a/4*hx;
    end
    for k = NR:N
        indx = N_c_c+(k-NR)*M;
        y_axis_s(i) = y_axis_s(i) + y(i, indx)*dr^3/6*bet_c*a_c*A_c/4*hx;
        for j = 2:M-1
            indx = N_c_c+(k-NR)*M-1+j;
            y_axis_s(i) = y_axis_s(i) + y(i, indx)/6*((r(j) + r(j+1))^3 - (r(j-1) + r(j))^3)*bet_c*a_c*A_c/4*hx;
        end
        indx = N_c_c+(k-NR)*M-1+M;
        y_axis_s(i) = y_axis_s(i) + y(i, indx)*4/3*(1 - (1 - dr/2)^3)*bet_c*a_c*A_c/4*hx;
    end
end
subplot(1, 3, 2);
plot(t, 100*(y_axis_s - y_axis_s(1))/y_axis_s(1), '-or', 'MarkerFaceColor', 'b');
%plot(t, y_axis_s, '-or', 'MarkerFaceColor', 'b');  % absolute values
xlim([0 t_max]);
title('Lithium Conservation in solid');
xlabel('Time [s]');
ylabel('Percentage deviation from initial value [%]');
grid on;

y_axis_sl = y_axis_liq + y_axis_s;
subplot(1, 3, 3);
plot(t, 100*(y_axis_sl - y_axis_sl(1))/y_axis_sl(1), '-or', 'MarkerFaceColor', 'r');
xlim([0 t_max]);
title('Total Lithium Conservation');
xlabel('Time [s]');
ylabel('Percentage deviation from initial value [%]');
grid on;
