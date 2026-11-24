function [q,V,res,iter,info] = morse_npss(q0,qmin,c)
%MORSE_NPSS Retained-history two-stage NPSS, fixed ritz-history scheme.
% Stage I works in the complement of a segment's spectral nullspace. A
% signed Rayleigh value below rayleigh_tol triggers a spectrum check;
% negative Rayleigh curvature switches to the full Hessian in Stage II.
% BB secants and transported three-block Ritz history are never reset.
q=q0(:); m=numel(q);
[~,~,Hmin]=morse_energy_gradient_hessian(qmin(:),c.a,c.r0);
[Nbar,Q]=spectral_basis(Hmin,c.null_tol);
assert(size(Nbar,2)==3,'Expected two translations and one rotation at the mother.');
[~,g,H]=morse_energy_gradient_hessian(q,c.a,c.r0);
V=lowest(H,Q); Vprevious=[];
info=struct('stage',1,'switch_iter',NaN,'switch_ritz',NaN,'null_updates',0, ...
    'null_checks',0,'bb_history_resets',0,'ritz_history_transfers',0,'status','maxiter');

% The same 1e-4 warm Euler step used by the other Table 1 methods.
q_old=q; D_old=g-2*V*(V'*g);
q=q-1e-4*D_old;
[~,g,H]=morse_energy_gradient_hessian(q,c.a,c.r0);
res=norm(g)/sqrt(m); bb_ready=true; iter=0;
while iter<c.maxiter
    if any(~isfinite(H(:))) || ~isfinite(res)
        info.status='nonfinite'; break;
    end
    if info.stage==1
        rho=V'*H*V;
        if rho<0
            % Release spectral constraints while transporting previous history.
            Qnext=eye(m);
            if ~isempty(Vprevious)
                Vprevious=Qnext'*(Q*Vprevious);
                info.ritz_history_transfers=info.ritz_history_transfers+1;
            end
            Q=Qnext; V=lowest(H,Q);
            info.stage=2; info.switch_iter=iter; info.switch_ritz=rho;
        end
        if info.stage==1
            if rho<c.rayleigh_tol
                info.null_checks=info.null_checks+1;
                [Nnew,Qnew]=spectral_basis(H,c.null_tol);
            else
                Nnew=Nbar;
            end
            if norm(Nnew*Nnew'-Nbar*Nbar','fro')>1e-6 && ~isempty(Qnew)
                if ~isempty(Vprevious)
                    Vprevious=Qnew'*(Q*Vprevious);
                    info.ritz_history_transfers=info.ritz_history_transfers+1;
                end
                Nbar=Nnew; Q=Qnew; V=lowest(H,Q);
                info.null_updates=info.null_updates+1;
            end
        end
    end
    if info.stage==2 && res<=c.res_tol
        info.status='converged'; break;
    elseif res==0
        info.status='stage1_stalled'; break;
    end

    % BB2 on the reflected force, preserving secants across all refreshes.
    D=g-2*V*(V'*g); alpha=c.dt;
    if bb_ready
        s=q-q_old; y=D-D_old; yy=y'*y; sy=s'*y;
        if yy>eps && isfinite(yy) && isfinite(sy), alpha=abs(sy/yy); end
    end
    alpha=min(max(alpha,c.bb_min_factor*c.dt),c.bb_max_factor*c.dt);
    q_old=q; D_old=D; bb_ready=true;
    q=q-alpha*D; iter=iter+1;
    [~,g,H]=morse_energy_gradient_hessian(q,c.a,c.r0);
    res=norm(g)/sqrt(m);

    % One Ritz step in the current allowed space, with transported history.
    if isfinite(res) && all(isfinite(H(:)))
        v=Q'*V;
        if norm(v)<=1e-12
            V=lowest(H,Q); continue;
        end
        A=Q'*H*Q;
        [v,Vprevious]=morse_ritz_step((A+A')/2,v,Vprevious);
        V=Q*v; V=V/norm(V);
    end
end
if ~isfinite(res) || any(~isfinite(H(:)))
    info.status='nonfinite';
elseif info.stage==2 && res<=c.res_tol
    info.status='converged';
end
end

function [N,Q]=spectral_basis(H,tol)
[Z,ev]=eig((H+H')/2,'vector');
isnull=abs(ev)<=tol; N=Z(:,isnull); Q=Z(:,~isnull);
assert(~isempty(Q),'The spectral nullspace leaves no ascent direction.');
end

function v=lowest(H,Q)
A=Q'*H*Q;
[Z,ev]=eig((A+A')/2,'vector'); [~,i]=min(ev);
v=Q*Z(:,i); v=v/norm(v);
end
