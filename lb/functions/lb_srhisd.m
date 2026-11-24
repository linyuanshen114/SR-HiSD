function [x,output]=lb_srhisd(F,H,x,options)
%LB_SRHISD FFT-cached SR-HiSD with the shared generalized-Gram LOBPCG step.
% One BB2/SI position step and one current-tangent Ritz step per iteration.
% F supports F(x,state_hat); H is the raw Hessian action. File I/O is owned
% by the caller via optional checkpoint_callback(state,history_table).
global LB
cfg=solver_settings(options); m=numel(x); k=cfg.k; clock=tic;
if isfield(options,'precd') && ~isempty(options.precd); precondition=options.precd;
else; precondition=@lb_precondition; end
C=empty_counts();
names={'iteration','force_rms','force_max','mean_phi','energy','bb_step','k', ...
 'ritz_min','direction_residual_max','orthogonality_error','null_leak_rms', ...
 'reference_rank','refresh_count','check_count','drift_curvature','force_count', ...
 'hvp_count','precondition_count','fft_forward','fft_inverse','tangent_count', ...
 'energy_count','wall_seconds','kernel_seconds','monitor_seconds','refresh_seconds', ...
 'io_seconds','peak_rss_kb','trigger_code','history_reset_count','bb_reset_count', ...
 'direction_seconds','direction_projection_seconds','direction_projection_calls','direction_projected_columns'};

%% Initialize the frame in the current translation-normal space.
rng(cfg.seed,'twister');
p=model('tangent',x,struct('tol',cfg.tangent_tol)); Q=p.Q;
V=lb_sr_orth(options.V(:,1:k),@(X) lb_sr_project(X,Q),1e-12);
assert(size(V,2)==k,'LB:InitialRank','Initial SR direction lost projected rank.');
S=struct('phi',x(:),'V',V,'Vprevious',zeros(m,0),'Nbar',zeros(m,0), ...
 'iteration',0,'xp',x(:),'gp',zeros(m,1),'bb_valid',false,'bta',cfg.dt, ...
 'theta',NaN(k,1),'kernel_seconds',0,'direction_seconds',0, ...
 'direction_projection_seconds',0,'direction_projection_calls',0, ...
 'direction_projected_columns',0,'io_seconds',0,'wall_seconds',0, ...
 'history',zeros(0,numel(names)));
S.model=struct('N',LB.N,'L',LB.L,'par',LB.par, ...
 'inner_product','X_transpose_Y_divided_by_N_cubed','state_constraint','mean_zero');
force_cache_valid=false; cached_force=[]; status='budget_exhausted';
history=zeros(cfg.maxiter,numel(names));
for it=1:cfg.maxiter
    tk=tic; x=S.phi; V=S.V;
    %% Cached force, reflected force and BB2 position step.
    if cfg.fft_cache_enabled && force_cache_valid
        f=cached_force; model('force_cache_hit');
    else
        f=model('force',x);
    end
    force_cache_valid=false;
    reflected=2*V*(V'*f/m); g=f-reflected; b=cfg.dt;
    if S.bb_valid
        dx=x-S.xp; dg=g-S.gp; denominator=dg'*dg/m;
        if abs(denominator)>1e-20; b=abs((dx'*dg/m)/denominator); end
        b=min(cfg.bb_max,max(cfg.bb_min,b));
    end
    S.xp=x; S.gp=g; S.bb_valid=true;
    x=model('position',x,reflected,b);

    %% Share one new-state FFT across tangent, force and energy operations.
    if cfg.fft_cache_enabled
        state_hat=model('state_spectrum',x);
        p=model('tangent_from_spectrum',x,struct('tol',cfg.tangent_tol),state_hat);
    else
        p=model('tangent',x,struct('tol',cfg.tangent_tol));
    end
    Q=p.Q; project=@(X) lb_sr_project(X,Q);
    direction_clock=tic;
    [V,Vprevious,info]=lb_sr_lobpcg_step(@(v) model('hv',x,v), ...
        V,S.Vprevious,project,@(v) model('precondition',v),cfg);
    S.direction_seconds=S.direction_seconds+toc(direction_clock);
    S.direction_projection_seconds=S.direction_projection_seconds+info.project_seconds;
    S.direction_projection_calls=S.direction_projection_calls+info.project_count;
    S.direction_projected_columns=S.direction_projected_columns+info.project_scalar_columns;
    S.phi=x; S.V=V; S.Vprevious=Vprevious; S.theta=info.ritz_values;
    S.bta=b; S.iteration=it;
    if cfg.fft_cache_enabled
        fnew=model('force_from_spectrum',x,state_hat);
        energy=model('energy_from_spectrum',x,state_hat);
        cached_force=fnew; force_cache_valid=true;
    else
        fnew=model('force',x); energy=model('energy',x);
    end
    residual=norm(fnew)/sqrt(m);
    S.kernel_seconds=S.kernel_seconds+toc(tk);
    assert(all(isfinite(x)) && isfinite(residual) && all(isfinite(V(:))), ...
        'LB:Nonfinite','Nonfinite SR state.');
    S.counts=C; S.rng=rng; S.wall_seconds=toc(clock);
    history(it,:)=[it residual max(abs(fnew)) mean(x) energy b k min(S.theta) ...
        max(info.residual_rms) norm(S.V'*S.V/m-eye(k),'fro') ...
        norm(Q'*S.V/m,'fro') size(Q,2) 0 0 NaN C.force C.hvp C.precondition ...
        C.fft_forward C.fft_inverse C.tangent C.energy S.wall_seconds ...
        S.kernel_seconds 0 0 S.io_seconds NaN 0 0 0 S.direction_seconds ...
        S.direction_projection_seconds S.direction_projection_calls S.direction_projected_columns];
    if mod(it,cfg.log_every)==0 || it==1 || residual<=cfg.tol
        fprintf('SR-HiSD it=%d RMS=%.6e theta=%+.6e rank=%d HVP=%d FFT=%d\n', ...
            it,residual,min(S.theta),size(Q,2),C.hvp,C.fft_forward+C.fft_inverse);
    end
    if isfield(options,'checkpoint_callback') && ~isempty(options.checkpoint_callback) && ...
            mod(it,cfg.checkpoint_every)==0
        S.history=history(1:it,:);
        io_clock=tic;
        options.checkpoint_callback(S,array2table(S.history,'VariableNames',names));
        S.io_seconds=S.io_seconds+toc(io_clock);
    end
    if residual<=cfg.tol; status='force_converged'; break; end
    if toc(clock)>=cfg.wall_seconds; status='wall_budget_exhausted'; break; end
end
S.history=history(1:S.iteration,:); S.rng=rng; S.wall_seconds=toc(clock);

%% Audit the endpoint with a fresh force and energy evaluation.
search_counts=C; audit_clock=tic;
force=model('force',S.phi); energy=model('energy',S.phi);
audit_seconds=toc(audit_clock); audit_counts=C;
for name=fieldnames(search_counts).'
    audit_counts.(name{1})=C.(name{1})-search_counts.(name{1});
end
S.counts=C; S.status=status; x=S.phi;
output=struct('x',x,'V',S.V,'it',S.iteration,'alph',S.theta, ...
    'res',norm(force)/sqrt(m),'force_max',max(abs(force)),'energy',energy, ...
    'status',status,'solve_seconds',S.kernel_seconds,'search_seconds',S.wall_seconds, ...
    'kernel_seconds',S.kernel_seconds,'direction_seconds',S.direction_seconds, ...
    'direction_projection_seconds',S.direction_projection_seconds, ...
    'direction_projection_calls',S.direction_projection_calls, ...
    'direction_projected_columns',S.direction_projected_columns,'io_seconds',S.io_seconds, ...
    'lobpcg_block_tol',cfg.block_tol,'fft_cache_enabled',cfg.fft_cache_enabled, ...
    'direction_kernel','shared_existing_projection_gram','counts',C, ...
    'search_counts',search_counts,'final_force_audit_counts',audit_counts, ...
    'final_force_audit_seconds',audit_seconds,'state',S,'settings',cfg, ...
    'history',array2table(S.history,'VariableNames',names));

    function varargout=model(action,varargin)
        % Count the same FFT/HVP operations as the staged cached solver.
        switch action
            case 'force'
                timer=tic; varargout{1}=F(varargin{1});
                C.force=C.force+1; C.fft_forward=C.fft_forward+1; C.fft_inverse=C.fft_inverse+1;
                C.seconds_force=C.seconds_force+toc(timer);
            case 'state_spectrum'
                timer=tic; phi=reshape(varargin{1},LB.N,LB.N,LB.N);
                phi=phi-mean(phi(:)); varargout{1}=fftn(phi);
                C.fft_forward=C.fft_forward+1; C.shared_state_fft=C.shared_state_fft+1;
                C.seconds_shared_state_fft=C.seconds_shared_state_fft+toc(timer);
            case 'force_from_spectrum'
                timer=tic; varargout{1}=F(varargin{1},varargin{2});
                C.force=C.force+1; C.fft_inverse=C.fft_inverse+1;
                C.seconds_force=C.seconds_force+toc(timer);
            case 'force_cache_hit'
                C.force_cache_hits=C.force_cache_hits+1;
            case 'hv'
                assert(size(varargin{2},2)==1,'HVP requires a single column.');
                timer=tic; varargout{1}=H(varargin{1},varargin{2});
                C.hvp=C.hvp+1; C.fft_forward=C.fft_forward+1; C.fft_inverse=C.fft_inverse+1;
                C.seconds_hvp=C.seconds_hvp+toc(timer);
            case 'precondition'
                assert(size(varargin{1},2)==1,'Preconditioner requires a single column.');
                timer=tic; varargout{1}=precondition(varargin{1});
                C.precondition=C.precondition+1; C.fft_forward=C.fft_forward+1; C.fft_inverse=C.fft_inverse+1;
                C.seconds_precondition=C.seconds_precondition+toc(timer);
            case 'energy'
                timer=tic; varargout{1}=lb_energy(varargin{1});
                C.energy=C.energy+1; C.fft_forward=C.fft_forward+1; C.fft_inverse=C.fft_inverse+1;
                C.seconds_energy=C.seconds_energy+toc(timer);
            case 'energy_from_spectrum'
                timer=tic; varargout{1}=lb_energy(varargin{1},varargin{2});
                C.energy=C.energy+1; C.fft_inverse=C.fft_inverse+1;
                C.seconds_energy=C.seconds_energy+toc(timer);
            case 'tangent'
                timer=tic; varargout{1}=lb_translation_projector(varargin{1},varargin{2});
                C.tangent=C.tangent+1; C.fft_forward=C.fft_forward+1; C.fft_inverse=C.fft_inverse+3;
                C.seconds_tangent=C.seconds_tangent+toc(timer);
            case 'tangent_from_spectrum'
                timer=tic; varargout{1}=lb_translation_projector(varargin{1},varargin{2},varargin{3});
                C.tangent=C.tangent+1; C.fft_inverse=C.fft_inverse+3;
                C.seconds_tangent=C.seconds_tangent+toc(timer);
            case 'position'
                timer=tic; phi=reshape(varargin{1},LB.N,LB.N,LB.N); phi=phi-mean(phi(:));
                reflected3=reshape(varargin{2},size(phi)); step=varargin{3}; par=LB.par;
                low=par.tau*phi-(par.gamma/2)*phi.^2+phi.^3/6;
                next=real(ifftn(fftn(phi+step*(-low-reflected3))./(LB.ones_hat+step*LB.Lhigh_hat)));
                next=next-mean(next(:)); varargout{1}=next(:);
                C.position=C.position+1; C.fft_forward=C.fft_forward+1; C.fft_inverse=C.fft_inverse+1;
                C.seconds_position=C.seconds_position+toc(timer);
        end
    end
end

function cfg=solver_settings(options)
cfg=struct();
defaults=struct('k',1,'maxiter',80000,'dt',.24,'tol',1e-7, ...
    'block_tol',1e-8,'gram_tol',100*eps,'tangent_tol',1e-6,'seed',20261028, ...
    'wall_seconds',43200,'log_every',100,'checkpoint_every',1000,'fft_cache_enabled',true);
for name=fieldnames(defaults).'
    if isfield(options,name{1}); cfg.(name{1})=options.(name{1});
    else; cfg.(name{1})=defaults.(name{1}); end
end
if ~isfield(options,'tol') && isfield(options,'epsf'); cfg.tol=options.epsf; end
if ~isfield(options,'log_every') && isfield(options,'outputp'); cfg.log_every=options.outputp; end
if ~isfield(options,'tangent_tol') && isfield(options,'proj_cfg') && isfield(options.proj_cfg,'tol')
    cfg.tangent_tol=options.proj_cfg.tol;
end
if isfield(options,'bb_min'); cfg.bb_min=options.bb_min;
elseif isfield(options,'betau'); cfg.bb_min=cfg.dt*options.betau;
else; cfg.bb_min=.2; end
if isfield(options,'bb_max'); cfg.bb_max=options.bb_max;
elseif isfield(options,'betat'); cfg.bb_max=cfg.dt*options.betat;
else; cfg.bb_max=2; end
assert(cfg.k>=1 && cfg.k==fix(cfg.k) && cfg.maxiter>=1);
assert(cfg.dt>0 && cfg.bb_min>0 && cfg.bb_max>=cfg.bb_min && cfg.log_every>=1);
assert(~(isfield(options,'nesterov') && options.nesterov) && ...
    ~(isfield(options,'heavy_ball') && options.heavy_ball),'LB:SRMomentum','Cached SR uses no momentum.');
if isfield(options,'direction_kernel')
    assert(strcmp(options.direction_kernel,'shared_existing_projection_gram'), ...
        'LB:SRDirectionKernel','Use the staged shared generalized-Gram direction kernel.');
end
end

function C=empty_counts()
C=struct('force',0,'hvp',0,'precondition',0,'energy',0,'tangent',0,'position',0, ...
    'fft_forward',0,'fft_inverse',0,'seconds_force',0,'seconds_hvp',0, ...
    'seconds_precondition',0,'seconds_energy',0,'seconds_tangent',0,'seconds_position',0, ...
    'shared_state_fft',0,'seconds_shared_state_fft',0,'force_cache_hits',0);
end
