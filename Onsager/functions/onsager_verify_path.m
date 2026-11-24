function validation = onsager_verify_path(route,p,cfg)
%ONSAGER_VERIFY_PATH Check actual residuals, constraints, lineage and the paper.
% No saved numerical reference result is an input to this verification.
assert(route.completed && all(diff([route.nodes.energy])<0));
mass_error = 0; seed_error = 0; residual_error = 0; frame_error = 0;
raw_error = 0; energy_error = 0;
for j = 1:numel(route.nodes)
    n = route.nodes(j);
    mass_error = max(mass_error,abs(sum(n.rho)*p.w-1));
    U = onsager_constraint_basis(n.logrho,p);
    [f,raw,energy] = onsager_force(n.logrho,U,p);
    residual_error = max(residual_error,abs(onsager_wnorm(f,p.w)-n.nf_projected));
    raw_error = max(raw_error,abs(onsager_wnorm(raw,p.w)-n.nf_raw));
    energy_error = max(energy_error,abs(energy-n.energy));
    frame_error = max(frame_error,norm(n.Vnegative'*n.Vnegative*p.w-eye(n.index),'fro'));
    assert(onsager_wnorm(f,p.w)<=cfg.tol && sum(n.eigenvalues < -cfg.index_tol)==n.index);
    if cfg.require_raw, assert(onsager_wnorm(raw,p.w)<=cfg.tol); end
    if j==1, continue; end
    q = find(arrayfun(@(a) a.edge.target_index==n.index,route.attempts),1,'last');
    a = route.attempts(q);
    parent = route.nodes(j-1);
    rho = max(parent.rho+a.solver.delta_rho,min(parent.rho));
    expected = log(rho/(sum(rho)*p.w));
    seed_error = max(seed_error,max(abs(expected-a.solver.initial_logrho)));
end
assert(mass_error<1e-12 && seed_error<1e-12 && residual_error<1e-11 && frame_error<1e-9);
assert(raw_error<1e-11 && energy_error<1e-10);
validation = struct('passed',true,'maximum_mass_error',mass_error, ...
    'maximum_negative_frame_orthogonality_error',frame_error, ...
    'maximum_seed_lineage_error',seed_error,'maximum_residual_recompute_error',residual_error, ...
    'maximum_raw_residual_recompute_error',raw_error,'maximum_energy_recompute_error',energy_error);
validation.paper=onsager_check_paper_result(route);
end

