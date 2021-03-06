function Ueq = Ueq_of_cs_anode(cs, cmax)
% Equilibrium potential for anode
% Ueq has unit of volts
csn = cs/cmax;

% LiC6 from either 
% (i) Srivinasan V., Newman J., (2004), "Design and optimization of a natural
% graphite/iron phosphate Lithium-ion cell". J. Electrochem. Soc., 151(10), A1530--A1538.
% (ii) Thomas K.E., (2002), "Thermal Modeling of Batteries with Porous Insertion
% Electrodes". University of Berkeley, PhD Thesis.
%cs_exp = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 ...
%          0.95 1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.62 1.65 1.68 1.7 ...
%          1.72 1.75 1.76 1.77]./1.77;
%V_exp  = [0.03 0.05 0.08 0.09 0.095 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.11 0.12 0.13 ...
%          0.13 0.135 0.135 0.135 0.135 0.135 0.135 0.135 0.15 0.16 0.165 0.17 0.175 0.185 0.195 0.205 ...
%          0.24 0.285 0.35 0.45 1 1.1];
%Ueq = interp1(1 - cs_exp, V_exp, real(csn));

% Analytical fit of the points above
%Ueq = 0.87643*exp(-80.7229*csn) + 0.0717649*exp(-11.1751*(csn-0.0515755)) - ...
%      0.00752717*tanh(178.2*(csn-0.238613)) - 0.0157348*tanh(30.019*(csn-0.529635)) - ...
%      0.931195*exp(22.9667*(csn-1.11242)) + 0.123513;

% Graphite anode from Ecker, Käbitz, Laresgoiti et al.
% Analytical fit (WebPlotDigitizer + gnuplot)
a = 0.716502;
b = 369.028;
c = 0.12193;
d = 35.6478;
e = 0.0530947;
g = 0.0169644;
h = 27.1365;
i = 0.312832;
j = 0.0199313;
k = 28.5697;
m = 0.614221;
n = 0.931153;
o = 36.328;
p = 1.10743;
q = 0.140031;
r = 0.0189193;
s = 21.1967;
t = 0.196176;

x = csn;

Ueq = a*exp(-b*x) + c*exp(-d*(x-e)) - r*tanh(s*(x-t)) - g*tanh(h*(x-i)) - j*tanh(k*(x-m)) - n*exp(o*(x-p)) + q;

% Constant for testing purpose
%Ueq = 0.2;
