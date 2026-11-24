function spectrum = lb_raw_spectrum(phi,opts)
%LB_RAW_SPECTRUM Matrix-free mean-zero low spectrum at one fixed LB field.
% All trial columns, residuals, and history lie in the allowed subspace.
% No full-space P*H*P eigensolve, dense complement, or grid-square matrix.
% Finite-start Ritz convergence is numerical evidence, not a full-spectrum
% interval proof. raw_residual_rms always uses the UNMODIFIED original H.
if nargin<2; opts=struct(); end
opts=defaults(opts,struct('k',8,'maxiter',1000,'res_tol',1e-8, ...
    'seed',20261028,'constraint_basis',[],'warm_V',[], ...
    'block_tol',1e-12,'check_every',10,'time_budget',Inf,'verbose',false, ...
    'direction_kernel','reference_generalized_gram'));
assert(strcmp(opts.direction_kernel,'reference_generalized_gram'), ...
    'LB:SpectrumKernel','Use the Gram low-spectrum direction kernel.');
phi=phi(:); m=numel(phi); t=tic;
assert(all(isfinite(phi)) && opts.k>=1 && opts.k==floor(opts.k));
old_rng=rng; restore_rng=onCleanup(@()rng(old_rng)); %#ok<NASGU>
rng(opts.seed,'twister');
Q=orth_columns(opts.constraint_basis,m,1e-11);
assert(opts.k < m-size(Q,2),'lb_raw_spectrum:Dimension','Too many requested eigenvectors.');
project=@(X)project_columns(X,Q,m);
V=orth_columns(project(opts.warm_V),m,opts.block_tol);
if size(V,2)>opts.k; V=V(:,1:opts.k); end
while size(V,2)<opts.k
    V=orth_columns([V,project(randn(m,opts.k-size(V,2)))],m,opts.block_tol);
end
Vp=[]; history=zeros(0,5); status='budget_exhausted';
last_info=struct(); iterations=0; hvps=0; preconditions=0;
try
    for iteration=1:opts.maxiter
        step=@lb_spectrum_lobpcg_step;
        [V,Vp,info]=step(@(v)lb_hessian_action(phi,v),V,Vp, ...
            project,@(v)lb_precondition(v),struct('block_tol',opts.block_tol));
        last_info=info; iterations=iteration;
        hvps=hvps+info.hvp_count; preconditions=preconditions+info.precondition_count;
        if iteration==1 || mod(iteration,opts.check_every)==0 || ...
                max(info.residual_rms)<=opts.res_tol || iteration==opts.maxiter
            history(end+1,:)=[iteration,toc(t),min(info.ritz_values), ...
                max(info.residual_rms),info.orthogonality_error]; %#ok<AGROW>
            if opts.verbose
                fprintf('spectrum it=%d residual=%.3e lowest=%+.10e\n', ...
                    iteration,max(info.residual_rms),min(info.ritz_values));
            end
        end
        if max(info.residual_rms)<=opts.res_tol
            status='projected_residual_converged'; break;
        end
        if toc(t)>=opts.time_budget; status='time_budget_exhausted'; break; end
    end
catch exception
    status=['failed:' exception.identifier];
    last_error=exception.message;
end
if ~exist('last_error','var'); last_error=''; end
% Independent fresh raw H applications, rather than cached search images.
HV=zeros(size(V));
for j=1:size(V,2); HV(:,j)=lb_hessian_action(phi,V(:,j)); end
hvps=hvps+size(V,2);
T=(V'*HV/m); T=(T+T')/2;
[rotation,lambda]=eig(T,'vector'); [lambda,order]=sort(real(lambda),'ascend');
rotation=rotation(:,order); V=V*rotation; HV=HV*rotation;
raw=HV-V*diag(lambda); constrained=project(HV)-V*diag(lambda);
raw_res=sqrt(sum(raw.^2,1)/m).'; proj_res=sqrt(sum(constrained.^2,1)/m).';
spectrum=struct('V',V,'lambda',lambda(:), ...
    'raw_residual_rms',raw_res,'projected_residual_rms',proj_res, ...
    'sign_interval_lower',lambda(:)-raw_res, ...
    'sign_interval_upper',lambda(:)+raw_res, ...
    'converged',all(proj_res<=opts.res_tol), ...
    'raw_converged',all(raw_res<=opts.res_tol), ...
    'status',status,'error',last_error,'iterations',iterations, ...
    'hvp_count',hvps,'precondition_count',preconditions,'seconds',toc(t), ...
    'constraint_rank',size(Q,2),'constant_explicitly_excluded',true, ...
    'mean_error',max(abs(mean(V,1))), ...
    'orthogonality_error',norm(V'*V/m-eye(size(V,2)),'fro'), ...
    'constraint_leak_rms',sqrt(sum((project(V)-V).^2,1)/m).', ...
    'constraint_overlap',Q'*V/m,'history',history,'options',opts, ...
    'last_direction_info',last_info, ...
    'completeness_note','Finite-start low-spectrum Ritz evidence; no rigorous full-spectrum no-missed-mode proof.');
end

function X=project_columns(X,Q,m)
if isempty(X); X=zeros(m,0); return; end
assert(size(X,1)==m,'lb_raw_spectrum:ConstraintSize','Wrong vector size.');
X=X-mean(X,1);
if ~isempty(Q); X=X-Q*(Q'*X/m); end
X=X-mean(X,1);
end

function Q=orth_columns(X,m,tol)
Q=zeros(m,0);
if isempty(X); return; end
assert(size(X,1)==m,'lb_raw_spectrum:BasisSize','Wrong basis size.');
for j=1:size(X,2)
    v=X(:,j)-mean(X(:,j));
    for pass=1:2; v=v-Q*(Q'*v/m); end
    n=norm(v)/sqrt(m);
    if n>tol; Q(:,end+1)=v/n; end %#ok<AGROW>
end
end

function a=defaults(a,b)
names=fieldnames(b);
for j=1:numel(names)
    if ~isfield(a,names{j}) || isempty(a.(names{j})); a.(names{j})=b.(names{j}); end
end
end
