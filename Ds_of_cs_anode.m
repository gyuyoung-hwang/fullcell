function Ds = Ds_of_cs_anode(cs, cmax)
% Diffusion coefficient in solid for anode, m^2/s

x = cs/cmax;

% Graphite anode from Ecker, Käbitz, Laresgoiti et al.
% Analytical fit (WebPlotDigitizer + gnuplot)

Dlog = 29021.6*x.^10 - 97916.1*x.^9 + 89874.2*x.^8 + 59572*x.^7 - 183309*x.^6 + 155891*x.^5 - ...
       66681*x.^4 + 15205.8*x.^3 - 1736.76*x.^2 + 78.0788*x - 9.85839;

Ds = (10.^Dlog)*1e-4;

% Constant for testing purpose
%Ds = 1e-14;
