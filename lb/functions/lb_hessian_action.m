function Hv = lb_hessian_action(phi, v)
% lb_hessian_action(phi, v)
% Hessian-vector product: D^2E(phi)[v] (projected to zero-mean)

global LB
N = LB.N;
par = LB.par;

phi3 = reshape(phi, N,N,N);
phi3 = phi3 - mean(phi3(:));

v3 = reshape(v, N,N,N);
v3 = v3 - mean(v3(:));

vhat = fftn(v3);

% High-order linear: xi^2*(Δ+q0^2)^2 v
Hv_high = real(ifftn(LB.Lhigh_hat .* vhat));

% Low-order/nonlinear:
% derivative of (tau*phi - gamma/2*phi^2 + 1/6*phi^3) wrt phi, applied to v:
% tau*v - gamma*phi*v + 1/2*phi^2*v
Hv_low = par.tau*v3 - par.gamma*(phi3.*v3) + 0.5*(phi3.^2).*v3;

Hv3 = Hv_high + Hv_low;
Hv3 = Hv3 - mean(Hv3(:));
Hv = Hv3(:);
end
