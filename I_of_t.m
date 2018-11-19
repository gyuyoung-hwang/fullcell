function I = I_of_t(t, I0, t0)
% Time dependent current draw
% t: current time
% t0: charge/discharge period, same units as t
% I0: maximum current draw, A

if t0 > 0
    % Sine function
    I = I0*sin(2*pi*t/t0);
    
    % Step function
    %if mod(floor(t/t0), 2) == 0  % even
    %    I = I0;
    %else  % odd
    %    I = -I0;
    %end
else
    I = I0;
end
