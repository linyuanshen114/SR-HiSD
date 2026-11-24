function node = onsager_analyze(logrho,target,p,cfg,iterations)
U = onsager_constraint_basis(logrho,p);
[f,raw,energy] = onsager_force(logrho,U,p);
Q = null(U')/sqrt(p.w);
HQ = onsager_hessian_action(logrho,Q,U,p);
G = Q'*HQ*p.w;
[eta,D] = eig((G+G')/2,'vector');
[D,order] = sort(real(D));
k = sum(D < -cfg.index_tol);
V = Q*real(eta(:,order(1:k)));
node = struct('id',sprintf('onsager_idx%d',target),'logrho',logrho, ...
    'rho',exp(logrho),'target_index',target,'index',k,'energy',energy, ...
    'nf_projected',onsager_wnorm(f,p.w),'nf_raw',onsager_wnorm(raw,p.w), ...
    'eigenvalues',D,'Vnegative',onsager_canonical_signs(V),'iterations',iterations, ...
    'accepted',k==target && onsager_wnorm(f,p.w)<=cfg.tol ...
        && (~cfg.require_raw || onsager_wnorm(raw,p.w)<=cfg.tol));
fprintf('CHECK index %d / %d, nf %.6e raw %.6e E %.12e accepted %d\n', ...
    k,target,node.nf_projected,node.nf_raw,energy,node.accepted);
end

