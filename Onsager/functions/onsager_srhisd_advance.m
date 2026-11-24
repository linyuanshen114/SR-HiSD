function s = onsager_srhisd_advance(s,opt,p,cfg,count)
U = onsager_constraint_basis(s.logrho,p);
[f,raw,~] = onsager_force(s.logrho,U,p);
for j = 1:count
    rho = exp(s.logrho);
    g = onsager_project(f-2*s.V*(s.V'*f*p.w),U,p.w);
    beta = opt.dt;
    if strcmp(opt.method,'bb')
        if s.it>0
            dg = g-s.gp;
            denom = dg'*dg*p.w;
            if denom > realmin
                beta = abs((s.Dx'*dg*p.w)/denom);
            end
        end
        beta = min(opt.beta_max,max(opt.beta_min,beta));
    end
    position_g = g;
    if ~cfg.project_position, position_g = g+raw-f; end
    rhs = rho+beta*s.logrho+beta*position_g;
    [next,implicit_error] = onsager_entropy_solve(beta,rhs,s.logrho);
    a = max(next);
    next = next-(a+log(sum(exp(next-a)))+log(p.w));
    nextU = onsager_constraint_basis(next,p);
    oldV = onsager_orthogonalize(onsager_project(s.V,nextU,p.w),p.w);
    if opt.target>0
        [s.V,s.Vp] = onsager_lobpcg_frame(next,oldV,s.Vp,nextU,p,opt.target);
    end
    s.gp = g;
    s.Dx = onsager_project(exp(next)-rho,nextU,p.w);
    s.logrho = next; U = nextU;
    [f,raw,s.energy] = onsager_force(next,U,p);
    s.it = s.it+1; s.nf = onsager_wnorm(f,p.w); s.raw_nf = onsager_wnorm(raw,p.w); s.beta = beta;
    s.history(end+1,:) = [s.it,s.nf,s.raw_nf,beta,s.energy,implicit_error,onsager_wnorm(s.Dx,p.w)];
    assert(all(isfinite(s.history(end,:))) && all(isfinite(next)),'Nonfinite iterate. Reduce the stage step.');
    if s.nf <= cfg.tol && (~cfg.require_raw || s.raw_nf <= cfg.tol), break; end
end
end

