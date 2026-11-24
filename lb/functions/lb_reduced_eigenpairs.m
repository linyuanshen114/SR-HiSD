function [V, lambda] = lb_reduced_eigenpairs(phi, k)
% lb_full_eigenpairs(phi,k)
% Use LOBPCG-style iteration (lb_full_ritz) to compute k smallest eigenpairs
% of Hessian at phi.

global LB
dof = LB.dof;

opt = lb_solver_options(1);
opt.k = k;
opt.maxiter = 100;
opt.precd = @lb_precondition;
opt.innpfunc = @(X,Y) (X'*Y) / dof;

rng(0);
opt.V = randn(dof,k);

[V, lambda] = lb_reduced_ritz(@lb_hessian_action, phi, opt);

% Orthonormalize
for i=1:k
    V(:,i) = V(:,i) - V(:,1:i-1)*(opt.innpfunc(V(:,1:i-1),V(:,i)));
    V(:,i) = V(:,i) / sqrt(opt.innpfunc(V(:,i),V(:,i)));
end

lambda = lambda(:);

end
