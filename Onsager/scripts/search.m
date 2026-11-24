%% Onsager SR-HiSD: follow the index 5 -> 4 -> 2 -> 0 path.
% Run this script directly. Full numerical settings are written below.
% To test two actual steps: review_options=struct('Quick',true); then Run.
task_dir=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(task_dir),'common'));
if ~exist('review_options','var'); review_options=struct(); end
ctx=review_setup(mfilename('fullpath'),'search',review_options);

%% 1. Define the discrete model and the three descent stages.
% E(rho)=w*sum(rho.*log(rho)) + (w/2)*rho'*(K*w)*rho, w*sum(rho)=1.
% K(p,q)=80*sqrt(1-(p.q)^2), with the inherited constant centering.
cfg=struct('data_dir',fullfile(ctx.data_dir,'input'), ...
    'tol',1e-7,'index_tol',1e-7,'lmax',12,'epsrot',1e-5,'sr_tol',1e-6, ...
    'project_position',false,'require_raw',true,'rng_seed',20261023);
stages=struct('target',{4,2,0},'method',{'constant','constant','constant'}, ...
    'dt',{0.03,0.03,0.5},'beta_min',{0.0002,0.001,0.01}, ...
    'beta_max',{0.03,0.03,2},'maxiter',{200000,120000,30000}, ...
    'modes',{5,[3 4],[1 2]},'signs',{-1,[-1 1],[1 1]}, ...
    'epsilon',{1e-4,1e-4,1e-3});
checkpoint_steps=1000;
if ctx.quick
    for level=1:3; stages(level).maxiter=2; end
    checkpoint_steps=1;
end
grid=load(fullfile(cfg.data_dir,'grid.mat'),'xyz');
seed=load(fullfile(cfg.data_dir,'seed.mat'),'logrho');
assert(isequal(size(grid.xyz),[2562 3]) && numel(seed.logrho)==2562);
assert(all(isfinite(grid.xyz(:))) && all(isfinite(seed.logrho(:))));
if strcmp(ctx.mode,'check')
    saved=load(fullfile(ctx.data_dir,'reference','result.mat'),'route');
    onsager_check_paper_result(saved.route);
    review_result=review_finish(ctx,struct('status','PASS','grid_points',2562));
    return;
end

%% 2. Assemble the kernel, quadrature and rotational constraint operators.
% constraint_basis removes mass and the three rotational orbit directions.
% hessian_action applies the quotient Hessian; no dense model is hidden in I/O.
rng(cfg.rng_seed,'twister');
p=onsager_make_problem(cfg);
active=[];
validation=struct();
if strcmp(ctx.mode,'reference')
    % Replay the current saved result, promoted from the validated full run.
    saved=load(fullfile(ctx.data_dir,'reference','result.mat'),'route');
    route=saved.route;
    validation=onsager_verify_path(route,p,cfg);
else
    %% 3. Analyze the seed, then perturb each accepted parent in released modes.
    first=onsager_analyze(seed.logrho,5,p,cfg,0);
    assert(first.accepted,'The supplied index-five seed failed certification.');
    route=struct('name','onsager','parameters',[0 0 80],'schedule',[5 4 2 0], ...
        'nodes',first,'edges',[],'attempts',[],'completed',false);
    for level=1:3
        opt=stages(level);
        parent=route.nodes(end);
        negative_frame=onsager_canonical_signs(parent.Vnegative);
        delta=opt.epsilon*negative_frame(:,opt.modes)*opt.signs(:);
        rho=max(parent.rho+delta,min(parent.rho));
        logrho=log(rho/(sum(rho)*p.w));
        active=onsager_initialize(logrho,negative_frame(:,1:opt.target),level,opt,p);
        active.parent_id=parent.id;
        active.delta_rho=delta;
        active.epsilon=opt.epsilon;
        active.signs=opt.signs;
        active.release_modes=opt.modes;
        active.settings{end+1}=struct('start_iteration',active.it,'options',opt);

        %% 4. Integrate SR-HiSD and the LOBPCG direction frame.
        % onsager_srhisd_advance explicitly calls force, the constraint
        % projector, the scalar entropy update and onsager_lobpcg_frame.
        fprintf('\nStage %d -> %d: dt=%g, budget=%d\n',parent.index,opt.target,opt.dt,opt.maxiter);
        while active.it<opt.maxiter && ...
                (active.nf>cfg.tol || (cfg.require_raw && active.raw_nf>cfg.tol))
            count=min(checkpoint_steps,opt.maxiter-active.it);
            active=onsager_srhisd_advance(active,opt,p,cfg,count);
            active.rng=rng;
            fprintf('%d: projected=%g, raw=%g, energy=%.12g\n', ...
                active.it,active.nf,active.raw_nf,active.energy);
            save(fullfile(ctx.run_dir,'checkpoint.mat'),'route','active','cfg','stages');
        end

        %% 5. Independently certify the endpoint before making it a parent.
        child=onsager_analyze(active.logrho,opt.target,p,cfg,active.it);
        edge=struct('from',parent.id,'to',child.id,'parent_index',parent.index, ...
            'target_index',opt.target,'child_index',child.index, ...
            'release_modes',active.release_modes,'signs',active.signs, ...
            'epsilon',active.epsilon,'method',opt.method,'dt',opt.dt, ...
            'iterations',active.it,'accepted',child.accepted,'selected',child.accepted);
        attempt=struct('node',child,'edge',edge,'solver',active);
        if isempty(route.attempts); route.attempts=attempt;
        else; route.attempts(end+1)=attempt; end
        history=array2table(active.history,'VariableNames', ...
            {'iteration','nf_projected','nf_raw','beta','energy','implicit_error','step_norm'});
        writetable(history,fullfile(ctx.run_dir,sprintf('history_index%d.csv',opt.target)));
        if ~child.accepted
            fprintf('Endpoint not accepted; the incomplete route and checkpoint are retained.\n');
            break;
        end
        route.nodes(end+1)=child;
        if isempty(route.edges); route.edges=edge; else; route.edges(end+1)=edge; end
        active=[];
    end
    route.completed=isequal([route.nodes.index],route.schedule) && all([route.nodes.accepted]);
    if route.completed; validation=onsager_verify_path(route,p,cfg); end
end

%% 6. Save the numerical path; plotting is a separate runnable script.
config=cfg;
save(fullfile(ctx.run_dir,'result.mat'),'route','active','config','stages','validation','-v7');
nodes=route.nodes;
node_table=table([nodes.index]',[nodes.energy]',[nodes.nf_projected]',[nodes.nf_raw]', ...
    'VariableNames',{'index','energy','projected_residual','raw_residual'});
writetable(node_table,fullfile(ctx.run_dir,'nodes.csv'));
status='PASS'; if ~route.completed; status='INCOMPLETE'; end
review_result=review_finish(ctx,struct('status',status,'completed',route.completed, ...
    'indices',[nodes.index],'numerical_search_performed',strcmp(ctx.mode,'compute')));
