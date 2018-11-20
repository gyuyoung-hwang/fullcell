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

% Constant for testing purpose
Ueq = 0.2;
