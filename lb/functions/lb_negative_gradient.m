function F = lb_negative_gradient(phi, phihat)
% lb_negative_gradient(phi)
% Returns negative gradient: F = - dE/dphi (projected to zero-mean subspace)

global LB
N = LB.N;
par = LB.par;

phi3 = reshape(phi, N,N,N);
phi3 = phi3 - mean(phi3(:));

if nargin < 2 || isempty(phihat); phihat = fftn(phi3); end

% High-order linear part: xi^2*(Δ+q0^2)^2 phi
grad_high = real(ifftn(LB.Lhigh_hat .* phihat));

% Low-order + nonlinear: tau*phi - gamma/2*phi^2 + 1/6*phi^3
grad_low = par.tau*phi3 - (par.gamma/2)*phi3.^2 + (1/6)*phi3.^3;

grad = grad_high + grad_low;

F3 = -grad;
F3 = F3 - mean(F3(:));     % projection P
F = F3(:);
end
