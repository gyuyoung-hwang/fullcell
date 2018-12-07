function Ueq = Ueq_of_cs_cathode(cs, cmax)
% Equilibrium potential for cathode
% Ueq has unit of volts
csn = cs/cmax;

% LFP (Srinivasan and Newman, J. Electrochem. Soc., 151(10):A1517--A1529, 2004)
%Ueq = (3.114559 + 4.438792*atan(-71.7352*csn + 70.85337) - 4.240252*atan(-68.5605*csn + 67.730082));

% NMC (Z. Mao and M. Farkhondeh, J. Electrochem. Soc.: Multi-particle model for a commercial blended...)
%Ueq = 6.51176 - 8.*csn + 7.1086.*csn.^2 - 1.55.*csn.^3 - 0.459.*csn.^6 - ...
%      5.00034e8.*exp(135.089.*csn.^2 - 118.089);

% NMC thick electrodes for Li-ion batteries (fit plot themselves)
%Ueq = 7.9760 - 5.5419.*csn + 5.2824.*csn.^1.07 - 4.0446.*csn.^0.0766 - ...
%      1.0556e-4.*exp(124.7407.*csn - 114.2593);

% LiNiCo from Ecker, Käbitz, Laresgoiti et al.
% Analytical fit (WebPlotDigitizer + gnuplot)
a = -2.35211;
c = 0.0747061;
d = 31.886;
e = 0.0219921;
g = 0.640243;
h = 5.48623;
i = 0.439245;
j = 3.82383;
k = 4.12167;
m = 0.176187;
n = 0.0542123;
o = 18.2919;
p = 0.762272;
q = 4.23285;
r = -6.34984;
s = 2.66395;
t = 0.174352;

x = csn;

Ueq = a*x - c*tanh(d*(x-e)) - r*tanh(s*(x-t)) - g*tanh(h*(x-i)) - j*tanh(k*(x-m)) - n*tanh(o*(x-p)) + q;

% Constant for testing purpose
%Ueq = 4.0;
