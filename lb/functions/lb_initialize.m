function lb_initialize(N, L, par)
% lb_initialize(N,L,par)
% Setup global environment for 3D LB model with FFT pseudo-spectral discretization.

global LB
LB = struct();
LB.N = N;
LB.L = L;
LB.par = par;
LB.dof = N^3;

% Grid in physical space: r_{ijk} = (L/N)*(i,j,k), i=0..N-1
x = (0:N-1) * (L/N);
LB.x = x;

% Precompute mesh for initial guesses / plotting
[X,Y,Z] = ndgrid(x, x, x);
LB.X = X; LB.Y = Y; LB.Z = Z;

% Fourier wavenumbers: q = (2*pi/L) * k, k = 0..N/2-1, -N/2..-1
k = [0:N/2-1, -N/2:-1];
q = (2*pi/L) * k;

[KX,KY,KZ] = ndgrid(q,q,q);
K2 = KX.^2 + KY.^2 + KZ.^2;
LB.K2 = K2;

k = (2*pi/L) * [0:(N/2-1), -N/2:-1];
[KX, KY, KZ] = ndgrid(k, k, k);

LB.KX = KX;
LB.KY = KY;
LB.KZ = KZ;


% High-order linear operator symbol: xi^2 * (q0^2 - |q|^2)^2
LB.Lhigh_hat = (par.xi^2) * (par.q0^2 - K2).^2;

% Denominator template for semi-implicit: 1 + dt*Lhigh_hat (dt changes per step)
LB.ones_hat = ones(N,N,N);

end
