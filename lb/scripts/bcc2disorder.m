%% BCC to disorder: native LBSD and FFT-cached SR-HiSD, exact kick 0.02
script_dir=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(script_dir)),'common'));
if ~exist('review_options','var'), review_options=struct(); end
ctx=review_setup(mfilename('fullpath'),'bcc2disorder',review_options);
if strcmp(ctx.mode,'check')
    assert(isfile(fullfile(ctx.data_dir,'arrays.h5')));
    assert(isfile(fullfile(ctx.data_dir,'input','bcc2disorder_srhisd.mat')));
    assert(isfile(fullfile(ctx.data_dir,'input','bcc2disorder_srhisd.json')));
    review_result=review_finish(ctx,struct('status','PASS','numerical_search_performed',false));
    return;
end

%% Periodic LB model and the selected disorder parameters
N=128; L=32*sqrt(2)*pi; kick=.02;
par=struct('tau',.05,'gamma',.6,'xi',1,'q0',1);
bb=struct('dt',.1,'betau',.2,'betat',2,'maxiter',80000,'epsf',1e-7,'outputp',100);
down_kick=.005; index_tolerance=1e-5;
attempt_rows={}; timing_rows={};
if strcmp(ctx.mode,'reference')
    routes=lb_reference_routes(ctx.data_dir,'bcc2disorder',kick);
else
    %% LBSD uses its own frozen original-space frame and exact phi_init.
    if ctx.quick
        seed=lb_demo_initial('bcc2disorder','lbsd4_to1',kick);
    else
        prefix='initial/bcc2disorder/kick_0020/lbsd4_to1/';
        seed=lb_read_view(ctx.data_dir,[prefix 'checkpoint_common_bcc.mat'], ...
            {'N','L','par','bb','kick','down_kick','phi_bcc','V_bcc_orig','V_bcc_tr'});
        exact=lb_read_view(ctx.data_dir,[prefix 'checkpoint_initial_search_sign_+1.mat'],{'phi_init'});
        seed.phi_init=exact.phi_init;
        assert(seed.N==N && abs(seed.L-L)<1e-12 && seed.kick==kick);
        assert(seed.par.tau==.05 && seed.par.gamma==.6 && seed.bb.dt==.1);
    end
    N=seed.N; bb=seed.bb; lb_initialize(N,seed.L,seed.par);
    innp=@(X,Y) (X'*Y)/N^3;
    v_shared=lb_normalize_rms(seed.V_bcc_tr(:,1),innp);
    V_native=seed.V_bcc_orig(:,1:4);
    v4=lb_normalize_rms(V_native(:,4),innp);
    if innp(v4,v_shared)<0, V_native(:,4)=-V_native(:,4); end
    rng(0,'twister'); method_clock=tic; verification_seconds=0;
    signs=[1 -1]; if ctx.quick, signs=1; end
    candidates={};

    %% Four-direction ascent; native index 1, or index 2 followed by descent.
    for sign=signs
        phi_start=seed.phi_bcc+sign*kick*v_shared;
        if sign==1, phi_start=seed.phi_init; end
        V_start=V_native; V_start(:,4)=sign*V_start(:,4);
        options=lb_options(4,V_start,bb,N);
        attempt_clock=tic;
        [phi,solver]=lb_lbsd(@lb_negative_gradient,@lb_hessian_action,phi_start,options);
        check_clock=tic; candidate=lb_candidate(phi,solver,[],index_tolerance,bb.epsf);
        check_seconds=toc(check_clock); verification_seconds=verification_seconds+check_seconds;
        candidate.valid=all(isfinite(phi)) && candidate.residual<=1.05*bb.epsf && ...
            ismember(candidate.index,[1 2]) && candidate.spectra_finite;
        row=lb_trial_record('lbsd4_to1','up_4',4,sign,candidate,toc(attempt_clock),check_seconds);
        attempt_rows{end+1}=row; candidates{end+1}=candidate; %#ok<SAGROW>
        lb_write_stage(ctx.run_dir,sprintf('lbsd_up4_sign_%+d',sign),candidate,row);
        if candidate.valid, break; end
    end
    current=lb_choose_candidate(candidates,2); stages={current.phi}; stage_indices=current.index;
    if current.valid && current.index==2
        negative=find(current.lambda_full < -index_tolerance);
        assert(numel(negative)==2,'LB:DownwardFrame','Native downward frame requires two negative modes.');
        V_keep=lb_normalize_rms(current.V_full(:,negative(1)),innp);
        v_drop=lb_normalize_rms(current.V_full(:,negative(2)),innp);
        candidates={}; signs=[1 -1]; if ctx.quick, signs=1; end
        for sign=signs
            phi_start=current.phi+sign*down_kick*v_drop;
            options=lb_options(1,V_keep,bb,N);
            attempt_clock=tic;
            [phi,solver]=lb_lbsd(@lb_negative_gradient,@lb_hessian_action,phi_start,options);
            check_clock=tic; candidate=lb_candidate(phi,solver,1,index_tolerance,bb.epsf);
            check_seconds=toc(check_clock); verification_seconds=verification_seconds+check_seconds;
            candidate.valid=candidate.valid && candidate.spectra_finite;
            row=lb_trial_record('lbsd4_to1','down_2_to_1',1,sign,candidate,toc(attempt_clock),check_seconds);
            attempt_rows{end+1}=row; candidates{end+1}=candidate; %#ok<SAGROW>
            lb_write_stage(ctx.run_dir,sprintf('lbsd_down2to1_sign_%+d',sign),candidate,row);
            if candidate.valid, break; end
        end
        current=lb_choose_candidate(candidates,1);
        if current.valid, stages{end+1}=current.phi; stage_indices(end+1)=current.index; end
    end
    search_seconds=toc(method_clock);

    %% Recompute the final diagnostics and follow both native GF branches.
    endpoints=struct([]); v1=[];
    if (current.valid && current.index==1) || ctx.quick
        check_clock=tic;
        [V_final_full,lambda_final_full]=lb_full_eigenpairs(current.phi,10); %#ok<NASGU>
        [V_final_reduced,lambda_final_reduced]=lb_reduced_eigenpairs(current.phi,10); %#ok<NASGU>
        v1=lb_normalize_rms(V_final_reduced(:,1),innp);
        verification_seconds=verification_seconds+toc(check_clock);
    end
    gf_seconds=0;
    if ~isempty(v1)
        gf=struct('dt',1,'maxiter',40000,'tol',1e-9,'outputEvery',200, ...
            'store_history',false,'nesterov',false,'chunk_iter',2000);
        if ctx.quick, gf.maxiter=3; end
        drawing=struct('draw_fig',false,'verbose',true);
        for branch=1:2
            branch_signs=[1 -1]; sign=branch_signs(branch);
            [phi_end,flow]=lb_flow_with_checkpoints(current.phi+sign*1e-5*v1,gf,drawing, ...
                ctx.run_dir,sprintf('lbsd_gf_sign_%+d',sign),current.phi,v1,sign);
            gf_seconds=gf_seconds+flow.elapsed_seconds; check_clock=tic;
            endpoint=lb_branch_info(phi_end,seed.phi_bcc,current.phi,sign);
            verification_seconds=verification_seconds+toc(check_clock);
            endpoint.flow=flow;
            if isempty(endpoints), endpoints=endpoint; else, endpoints(branch)=endpoint; end
        end
    end
    routes=lb_make_route('bcc2disorder','lbsd4_to1',seed,current,stages,endpoints,ctx.quick,stage_indices);
    lbsd_rows=vertcat(attempt_rows{:});
    timing_rows{1}=struct('method',"lbsd4_to1",'version',"native_blocktol_1e8", ...
        'iterations',sum([lbsd_rows.iterations]),'search_seconds',search_seconds, ...
        'kernel_seconds',sum([lbsd_rows.kernel_seconds]),'verification_seconds',verification_seconds, ...
        'gf_seconds',gf_seconds,'elapsed_before_plot_seconds',toc(method_clock),'complete',routes.complete);

    %% Cached SR-HiSD 32410: exact shared input and current-tangent Ritz frame.
    if ctx.quick
        seed=lb_demo_initial('bcc2disorder','srhisd1',kick);
        V_seed=seed.V_bcc_tr;
    else
        stored=load(fullfile(ctx.data_dir,'input','bcc2disorder_srhisd.mat'),'input'); input=stored.input;
        assert(input.N==128 && abs(input.L-32*sqrt(2)*pi)<1e-12);
        assert(input.par.tau==.05 && input.par.gamma==.6);
        seed=struct('N',input.N,'L',input.L,'par',input.par,'kick',.02, ...
            'phi_bcc',input.phi_parent(:),'phi_init',input.phi_init(:), ...
            'bb',struct('dt',.1,'betau',.2,'betat',2,'maxiter',80000,'epsf',1e-7,'outputp',100));
        V_seed=input.V_seed(:,1);
    end
    N=seed.N; lb_initialize(N,seed.L,seed.par);
    options=lb_options(1,V_seed,seed.bb,N);
    options.seed=20261028; options.tangent_tol=1e-6; options.block_tol=1e-8;
    options.gram_tol=100*eps; options.direction_kernel='shared_existing_projection_gram';
    options.fft_cache_enabled=true; options.wall_seconds=43200;
    options.checkpoint_every=1000;
    options.checkpoint_callback=@(state,history) lb_write_solver_checkpoint(ctx.run_dir,'srhisd1',state,history);
    method_clock=tic;
    [phi,solver]=lb_srhisd(@lb_negative_gradient,@lb_hessian_action,seed.phi_init,options);
    writetable(solver.history,fullfile(ctx.run_dir,'srhisd1_history.csv'));
    lb_write_solver_checkpoint(ctx.run_dir,'srhisd1',solver.state,solver.history);

    %% Verify a fresh raw negative direction; no index-one claim from force alone.
    force_converged=strcmp(solver.status,'force_converged'); negative_verified=false;
    spectrum=struct(); check_clock=tic;
    if force_converged || ctx.quick
        spectrum_options=struct('k',1,'maxiter',2500,'res_tol',1e-8,'seed',20261028, ...
            'constraint_basis',[],'warm_V',solver.V,'block_tol',1e-10,'time_budget',7200,'verbose',true);
        if ctx.quick, spectrum_options.maxiter=4; spectrum_options.time_budget=30; end
        spectrum=lb_raw_spectrum(phi,spectrum_options);
        negative_verified=~startsWith(spectrum.status,'failed') && ...
            spectrum.raw_residual_rms(1)<=1e-8 && ...
            spectrum.lambda(1)+spectrum.raw_residual_rms(1)<-1e-5;
    end
    verification_seconds=toc(check_clock);
    current=struct('phi',phi,'out',solver,'residual',solver.res,'energy',solver.energy, ...
        'index',NaN,'valid',force_converged && negative_verified, ...
        'negative_mode_verified',negative_verified,'strict_index_certified',false);
    row=lb_trial_record('srhisd1','direct_sr',1,1,current,solver.search_seconds,verification_seconds);
    attempt_rows{end+1}=row;
    lb_write_stage(ctx.run_dir,'srhisd1_final',current,row);
    save(fullfile(ctx.run_dir,'srhisd1_negative_mode.mat'),'spectrum','verification_seconds','-v7.3');
    endpoints=struct([]); gf_seconds=0;
    if current.valid || ctx.quick
        if negative_verified, v1=spectrum.V(:,1); else, v1=solver.V(:,1); end
        innp=@(X,Y) (X'*Y)/N^3; v1=lb_normalize_rms(v1,innp);
        gf=struct('dt',1,'maxiter',40000,'tol',1e-9,'outputEvery',200, ...
            'store_history',false,'nesterov',false,'chunk_iter',2000);
        if ctx.quick, gf.maxiter=3; end
        drawing=struct('draw_fig',false,'verbose',true);
        for branch=1:2
            branch_signs=[1 -1]; sign=branch_signs(branch);
            [phi_end,flow]=lb_flow_with_checkpoints(phi+sign*1e-5*v1,gf,drawing, ...
                ctx.run_dir,sprintf('srhisd1_gf_sign_%+d',sign),phi,v1,sign);
            gf_seconds=gf_seconds+flow.elapsed_seconds; check_clock=tic;
            endpoint=lb_branch_info(phi_end,seed.phi_bcc,phi,sign);
            verification_seconds=verification_seconds+toc(check_clock);
            endpoint.flow=flow;
            if isempty(endpoints), endpoints=endpoint; else, endpoints(branch)=endpoint; end
        end
    end
    routes(2)=lb_make_route('bcc2disorder','srhisd1',seed,current,{phi},endpoints,ctx.quick,1);
    timing_rows{2}=struct('method',"srhisd1",'version',"sr32410_fft_cache", ...
        'iterations',solver.it,'search_seconds',solver.search_seconds,'kernel_seconds',solver.kernel_seconds, ...
        'verification_seconds',verification_seconds,'gf_seconds',gf_seconds, ...
        'elapsed_before_plot_seconds',toc(method_clock),'complete',routes(2).complete);
end

%% Save the selected four-case version, measurements and automatic case plots.
save(fullfile(ctx.run_dir,'result.mat'),'routes','-v7.3');
if ~isempty(attempt_rows)
    Trials=struct2table(vertcat(attempt_rows{:})); Timing=struct2table(vertcat(timing_rows{:}));
    writetable(Trials,fullfile(ctx.run_dir,'attempts.csv'));
    writetable(Timing,fullfile(ctx.run_dir,'timing.csv'));
end
completed=all([routes.complete]); status='PASS'; if ~completed, status='INCOMPLETE'; end
figure_files={}; figure_records=struct([]);
if any([routes.plot_ready])
    [figure_files,figure_records]=lb_plot_case_routes(routes,ctx.run_dir,ctx.visible);
    save(fullfile(ctx.run_dir,'result.mat'),'figure_files','figure_records','-append');
end
review_result=review_finish(ctx,struct('status',status,'completed',completed,'route_count',numel(routes), ...
    'kick',kick,'numerical_search_performed',strcmp(ctx.mode,'compute'), ...
    'parameter_version','20261120','sr_implementation','sr32410_fft_cache','figure_files',{figure_files}));
