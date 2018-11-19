function yt = scheme_FEM(t, y, x, r)
% Computational scheme - Finite Element Method

parameters

rx = r.*a;
I_t = I_of_t(t, I0, t0);

% Zero all
yt   = zeros(N_end, 1);     % Vector-function
G    = zeros(N, 1);         % Reaction rate * bet (cell-centered)
Gn   = zeros(N, 1);         % Reaction rate (nodal)
csa  = zeros(N, 1);         % Concentration in solid at the boundary: c_s(r = a)
phis = zeros(N, 1);         % Potential in solid: phi_s(x,r)
qc1  = zeros(N, 1);         % Helper vector for concentration in liquid equation
qc2  = zeros(N, 1);         % Helper vector for concentration in liquid equation
qf1  = zeros(N, 1);         % Helper vector for potential in liquid equation
qf2  = zeros(N, 1);         % Helper vector for potential in liquid equation

% Extracting parameters
c      = y(1:N);
phi    = y(N+1:2*N);
phis_a = y(2*N+1:2*N+NL);
phis_c = y(2*N+NL+1:3*N+1-Ndelta);

% Helper vectors (nodal)
for i = 1:NL  % Anode
    csa(i) = y(N_c_a + i*M - 1);
    phis(i) = phis_a(i);
    eta = phi(i) - phis(i) + Ueq_of_cs_anode(csa(i), cmax_a);
    value = F*eta/(2*R*T);
    Gn(i) = k0_a*sqrt(c(i)*csa(i)*(cmax_a - csa(i)))*(exp(-value) - exp(value));
end
for i = NR:N  % Cathode
    j = i - NR + 1;  % starts from 1
    csa(i) = y(N_c_c + j*M - 1);
    phis(i) = phis_c(j);
    eta = phi(i) - phis(i) + Ueq_of_cs_cathode(csa(i), cmax_c);
    value = F*eta/(2*R*T);
    Gn(i) = k0_c*sqrt(c(i)*csa(i)*(cmax_c - csa(i)))*(exp(-value) - exp(value));
end

% Helper vectors (cell-centered)
for i = 1:N-1
    c_mid = (c(i) + c(i+1))*0.5;
    kappa_mid = kappa_of_c(c_mid);
    rtf_term = 2*R*T*(1 - tplus)/(c_mid*F^2);

    qc1(i) = -B(i)*(D_of_c(c_mid) + kappa_mid*rtf_term*(1 - tplus));
    qc2(i) = B(i)*kappa_mid*(1 - tplus)/F;

    qf1(i) = -B(i)*kappa_mid/F;
    qf2(i) = B(i)*kappa_mid*rtf_term;

    csa_mid = (csa(i) + csa(i+1))*0.5;
    eta = (phi(i) + phi(i+1))*0.5 - (phis(i) + phis(i+1))*0.5;
    if i < NL  % Anode
        eta = eta + Ueq_of_cs_anode(csa_mid, cmax_a);
    elseif i >= NR  % Cathode
        eta = eta + Ueq_of_cs_cathode(csa_mid, cmax_c);
    else  % Separator
        eta = 0;
    end
    value = F*eta/(2*R*T);
    G(i) = bet(i)*k0(i)*sqrt(c_mid*csa_mid*(cmax(i) - csa_mid))*(exp(-value) - exp(value));  % G*bet
end

% Concentration in liquid
yt(1) = -(qc1(1)*(c(2) - c(1)) + qc2(1)*(phi(2) - phi(1)))/(el(1)*hx^2);
for i = 2:N-1
    yt(i) = (qc1(i-1)*(c(i) - c(i-1)) - qc1(i)*(c(i+1) - c(i)) + ...
             qc2(i-1)*(phi(i) - phi(i-1)) - qc2(i)*(phi(i+1) - phi(i)))/(el(i)*hx^2);
end
yt(N) = (qc1(N-1)*(c(N) - c(N-1)) + qc2(N-1)*(phi(N) - phi(N-1)))/(el(N)*hx^2);

% Potential in liquid
yt(N+1) = phi(1);  % Dirichlet BC
for i = 2:N-1
    yt(N+i) = (qf1(i-1)*(phi(i) - phi(i-1)) - qf1(i)*(phi(i+1) - phi(i)) + ...
               qf2(i-1)*(c(i) - c(i-1)) - qf2(i)*(c(i+1) - c(i)))/hx + hx*(G(i-1) + G(i))*0.5;
end
yt(2*N) = (qf1(N-1)*(phi(N) - phi(N-1)) + qf2(N-1)*(c(N) - c(N-1)))/hx + hx*G(N-1)*0.5;

% Potential in solid (Anode)
yt(2*N+1) = I_t/A + sigma_s_a*(phis(2) - phis(1))/hx - F*hx*G(1)*0.5;
for i = 2:NL-1
    yt(2*N+i) = sigma_s_a*(phis(i+1) - 2*phis(i) + phis(i-1))/hx - F*hx*(G(i-1) + G(i))*0.5;
end
yt(2*N+NL) = -sigma_s_a*(phis(NL) - phis(NL-1))/hx - F*hx*G(NL-1)*0.5;

% Potential in solid (Cathode)
yt(2*N+NL+1) = sigma_s_c*(phis(NR+1) - phis(NR))/hx - F*hx*G(NR)*0.5;
for i = NR+1:N-1
    j = i - NR + 1;  % starts from 2
    yt(2*N+NL+j) = sigma_s_c*(phis(i+1) - 2*phis(i) + phis(i-1))/hx - F*hx*(G(i-1) + G(i))*0.5;
end
yt(3*N+1-Ndelta) = -I_t/A - sigma_s_c*(phis(N) - phis(N-1))/hx - F*hx*G(N-1)*0.5;

% Concentration in solid (Anode)
for i = 1:NL
    indx = N_c_a+(i-1)*M;
    yt(indx) = 0.25*Ds_a*(y(indx+1) - y(indx))/hr^2;
    for j = 2:M-1
        indx = N_c_a+(i-1)*M-1+j;
        yt(indx) = Ds_a*0.25*((rx(j+1) + rx(j))^2*(y(indx+1) - y(indx)) - ...
                              (rx(j-1) + rx(j))^2*(y(indx) - y(indx-1)))/hr^4;
    end
    indx = N_c_a+i*M-1;
    yt(indx) = -a^2*Gn(i)/hr^3 - (a - 0.5*hr)^2/hr^4*Ds_a*(y(indx) - y(indx-1));
end

% Concentration in solid (Cathode)
for i = NR:N
    indx = N_c_c+(i-NR)*M;
    yt(indx) = 0.25*Ds_c*(y(indx+1) - y(indx))/hr^2;
    for j = 2:M-1
        indx = N_c_c+(i-NR)*M-1+j;
        yt(indx) = Ds_c*0.25*((rx(j+1) + rx(j))^2*(y(indx+1) - y(indx)) - ...
                              (rx(j-1) + rx(j))^2*(y(indx) - y(indx-1)))/hr^4;
    end
    indx = N_c_c+(i+1-NR)*M-1;
    yt(indx) = -a^2*Gn(i)/hr^3 - (a - 0.5*hr)^2/hr^4*Ds_c*(y(indx) - y(indx-1));
end
