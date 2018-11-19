function Ueq = Ueq_of_cs_cathode(cs, cmax)
% Equilibrium potential for cathode
% Ueq has unit of volts
csn = cs/cmax;

% LFP (Srinivasan and Newman, J. Electrochem. Soc., 151(10):A1517--A1529, 2004)
%Ueq = (3.114559 + 4.438792*atan(-71.7352*csn + 70.85337) - 4.240252*atan(-68.5605*csn + 67.730082));

% NMC (Z. Mao and M. Farkhondeh, J. Electrochem. Soc.: Multi-particle model for a commercial blended...)
Ueq = 6.51176 - 8.*csn + 7.1086.*csn.^2 - 1.55.*csn.^3 - 0.459.*csn.^6 - ...
      5.00034e8.*exp(135.089.*csn.^2 - 118.089);

% NMC thick electrodes for Li-ion batteries (fit plot themselves)
%Ueq = 7.9760 - 5.5419.*csn + 5.2824.*csn.^1.07 - 4.0446.*csn.^0.0766 - ...
%      1.0556e-4.*exp(124.7407.*csn - 114.2593);

% Constant for testing purpose
%Ueq = 4.0;
