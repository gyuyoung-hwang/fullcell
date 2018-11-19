function Mass = mass_matrix_FEM(t, y, N, M, Ndelta, N_end)
% Mass matrix for Finite Element Method

dr = 1/(M - 1);
r = linspace(0, 1, M);

MN = eye(N, N);
MM = eye(M, M);

MN(1,1) = 1/3;
for i = 2:N
    MN(i,i) = 2/3;
    MN(i-1,i) = 1/6;
    MN(i,i-1) = 1/6;
end
MN(N,N) = 1/3;  % override the last element

MM(1,1) = 1/30;
for i = 2:M
    MM(i,i) = 1/15 + 2/3*(r(i)/dr)^2;
    MM(i-1,i) = 1/20 + 1/6*r(i)*r(i-1)/dr^2;
    MM(i,i-1) = MM(i-1,i);
end
MM(M,M) = 1/30 + 1/(3*dr^2) - 1/(6*dr);  % override the last element

Mass = eye(N_end);
Mass(1:N, 1:N) = MN;
Mass(N+1:3*N+1-Ndelta, N+1:3*N+1-Ndelta) = 0;
for i = 1:N-Ndelta+1
    Mass(3*N+2-Ndelta+(i-1)*M:3*N+1-Ndelta+i*M, 3*N+2-Ndelta+(i-1)*M:3*N+1-Ndelta+i*M) = MM;
end
Mass = sparse(Mass);
