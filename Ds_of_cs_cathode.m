function Ds = Ds_of_cs_cathode(cs, cmax)
% Diffusion coefficient in solid for cathode, m^2/s

x = cs/cmax;

% Graphite anode from Ecker, Käbitz, Laresgoiti et al.
% Analytical fit (WebPlotDigitizer + gnuplot)

a = 2.18899;
b = 11.6408;
c = 0.277226;
d = 1.97314;
e = 57.0596;
g = 0.817451;
h = -10.5785;

Dlog = a*exp(-b*(x-c).*(x-c)) + d*exp(-e*(x-g).*(x-g)) + h;

Ds = (10.^Dlog)*1e-4;

% Constant for testing purpose
%Ds = 1e-14;
