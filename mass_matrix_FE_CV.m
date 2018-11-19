function Mass = mass_matrix_FE_CV(t, y, N, M, Ndelta, N_end)
% Mass matrix for Finite Element + Control Volume Method

dr = 1/(M - 1);
r = linspace(0, 1, M);

MN = eye(N, N);
M1 = eye(M, M);
M2 = eye(M, M);

MN(1,1) = 1/3;
for i = 2:N
    MN(i,i) = 2/3;
    MN(i-1,i) = 1/6;
    MN(i,i-1) = 1/6;
end
MN(N,N) = 1/3;  % override the last element

M1(1,1) = 4/3*pi*(dr/2)^3;
for i = 2:M-1
    M1(i,i) = 4*pi*(r(i)^2*dr + dr^3/12);
end
M1(M,M) = 4/3*pi*(1 - (1 - dr/2)^3);

M2(1,1) = 3/4;
M2(1,2) = 1/8;
M2(2,1) = 1/4;
M2(2,2) = 6/8;
M2(2,3) = 1/8;
for i = 3:M-2
    M2(i,i) = 6/8;
    M2(i,i-1) = 1/8;
    M2(i,i+1) = 1/8;
end
M2(M-1,M-2) = 1/8;
M2(M-1,M-1) = 6/8;
M2(M-1,M) = 1/4;
M2(M,M-1) = 1/8;
M2(M,M) = 3/4;

MM = M2*M1;

Mass = eye(N_end);
Mass(1:N, 1:N) = MN;
Mass(N+1:3*N+1-Ndelta, N+1:3*N+1-Ndelta) = 0;
for i = 1:N-Ndelta+1
    Mass(3*N+2-Ndelta+(i-1)*M:3*N+1-Ndelta+i*M, 3*N+2-Ndelta+(i-1)*M:3*N+1-Ndelta+i*M) = MM;
end
Mass = sparse(Mass);
