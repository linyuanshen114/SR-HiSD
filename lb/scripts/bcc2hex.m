%% BCC to HEX: historical-input sequential LBSD and FFT-cached SR-HiSD
script_dir=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(script_dir)),'common'));
if ~exist('review_options','var'), review_options=struct(); end
ctx=review_setup(mfilename('fullpath'),'bcc2hex',review_options);
if strcmp(ctx.mode,'check')
    assert(isfile(fullfile(ctx.data_dir,'arrays.h5')));
    assert(isfile(fullfile(ctx.data_dir,'input','bcc2hex_lbsd.mat')));
    assert(isfile(fullfile(ctx.data_dir,'input','bcc2hex_srhisd.mat')));
    assert(isfile(fullfile(ctx.data_dir,'input','bcc2hex_srhisd.json')));
    review_result=review_finish(ctx,struct('status','PASS','numerical_search_performed',false));
    return;
end

%% Model and the four-stage original-space LBSD parameters
N=128; L=24*sqrt(2)*pi;
par=struct('tau',-.3,'gamma',.8,'xi',1,'q0',1);
bb=struct('dt',.24,'betau',.2/.24,'betat',2/.24, ...
    'maxiter',80000,'epsf',1e-7,'outputp',100);
kick=.01; down_kick=.001; kick_schedule=[.001 .002 .004];
attempt_rows={}; timing_rows={};
if strcmp(ctx.mode,'reference')
    routes=lb_reference_routes(ctx.data_dir,'bcc2hex');
else
    hex_reference=[];
    if ~ctx.quick
        hex_reference=lb_read_view(ctx.data_dir, ...
            'reruns/new_hex_gf_20261024_231111_w2esf5aw/hex_phase_certification_20261025/reference_hex.mat', ...
            {'phi','N','L','par'});
        assert(hex_reference.N==128 && abs(hex_reference.L-L)<1e-12);
    end
    %% Use the saved historical solver entrance, including all four columns.
    if ctx.quick
        seed=lb_demo_initial('bcc2hex','lbsd4_down',kick);
        V_start=seed.V_bcc_orig(:,1:4);
    else
        H=load(fullfile(ctx.data_dir,'input','bcc2hex_lbsd.mat'), ...
            'phi_init','V_init','bb','initial_alignment');
        common=lb_read_view(ctx.data_dir,'initial/bcc2hex/seed_common/checkpoint_common_bcc.mat',{'phi_bcc'});
        assert(numel(H.phi_init)==N^3 && isequal(size(H.V_init),[N^3 4]));
        assert(H.bb.dt==.24 && H.bb.min_step==.2 && H.bb.max_step==2 && H.bb.epsf==1e-7);
        seed=struct('N',N,'L',L,'par',par,'bb',bb,'kick',kick,'down_kick',down_kick, ...
            'phi_bcc',common.phi_bcc,'phi_init',H.phi_init(:));
        V_start=H.V_init;  % Pass unchanged; the LBSD kernel performs its own GS.
        clear H common
    end
    N=seed.N; L=seed.L; par=seed.par; bb=seed.bb;
    lb_initialize(N,L,par); innp=@(X,Y) (X'*Y)/N^3;
    method_clock=tic; verification_seconds=0;

    %% Four-direction ascent and independent warm/cold original-Hessian gate.
    options=lb_options(4,V_start,bb,N);
    attempt_clock=tic;
    initial_force=lb_negative_gradient(seed.phi_init);
    initial_rms=norm(initial_force)/sqrt(N^3);
    [phi,solver]=lb_lbsd(@lb_negative_gradient,@lb_hessian_action,seed.phi_init,options);
    check_clock=tic;
    current=lb_original_candidate(phi,solver,4,ctx.run_dir,'up_4',ctx.quick);
    check_seconds=toc(check_clock); verification_seconds=verification_seconds+check_seconds;
    row=lb_trial_record('lbsd4_down','up_4',4,1,current,toc(attempt_clock),check_seconds);
    row.initial_rms=initial_rms; row.kick_rms=0; row.initial_drive=NaN;
    row.return_parent_rms=NaN; row.gate_status=string(current.gate.status);
    attempt_rows{end+1}=row;
    lb_write_stage(ctx.run_dir,'lbsd_up4',current,row);
    stages={current.phi}; stage_indices=current.index;

    %% Follow 4 -> 3 -> 2 -> 1; try both signs before changing kick size.
    % Increase the kick only if both branches were ineffective or returned
    % to the verified parent. An unresolved or exhausted branch is retained.
    for target=[3 2 1]
        if ~current.valid || current.index~=target+1, break; end
        source=current.index;
        ids=find(current.gate.signs(:)==-1);
        assert(numel(ids)==source,'LB:DownwardFrame','Negative raw-HVP prefix disagrees with its index.');
        assert(lb_next_lbsd_k(source)==target);
        drop=current.V_full(:,ids(end)); drop=drop-mean(drop);
        drop=drop/sqrt(innp(drop,drop));
        V_keep=current.V_full(:,ids(1:target));
        parent_force=lb_negative_gradient(current.phi);
        accepted=false;
        for step_kick=kick_schedule
            fallback=false(1,2); signs=[1 -1];
            for branch=1:2
                sign=signs(branch);
                phi_start=current.phi+sign*step_kick*drop;
                initial_force=lb_negative_gradient(phi_start);
                initial_rms=norm(initial_force)/sqrt(N^3);
                kick_rms=sqrt(innp(phi_start-current.phi,phi_start-current.phi));
                initial_drive=innp(drop,initial_force-parent_force);
                label=sprintf('down_%dto%d_kick_%04d_sign_%+d',source,target,round(step_kick*1e6),sign);
                fprintf('DOWN %d->%d sign %+d kick %.8g force %.8g driving %.8g\n', ...
                    source,target,sign,kick_rms,initial_rms,initial_drive);
                if initial_rms<=bb.epsf
                    fallback(branch)=lb_kick_can_increase(true,false,NaN,source,NaN,step_kick,NaN,current.energy);
                    % Record a skipped, ineffective perturbation without a solve.
                    skipped=current; skipped.phi=phi_start; skipped.valid=false;
                    skipped.residual=initial_rms; skipped.energy=lb_energy(phi_start); skipped.index=NaN;
                    skipped.out=struct('it',0,'V',V_keep,'solve_seconds',0);
                    row=lb_trial_record('lbsd4_down',label,target,sign,skipped,0,0);
                    row.initial_rms=initial_rms; row.kick_rms=kick_rms; row.initial_drive=initial_drive;
                    row.return_parent_rms=NaN; row.gate_status="kick_below_tolerance";
                    attempt_rows{end+1}=row; %#ok<SAGROW>
                    continue;
                end
                options=lb_options(target,V_keep,bb,N);
                attempt_clock=tic;
                [phi,solver]=lb_lbsd(@lb_negative_gradient,@lb_hessian_action,phi_start,options);
                check_clock=tic;
                candidate=lb_original_candidate(phi,solver,target,ctx.run_dir,label,ctx.quick);
                check_seconds=toc(check_clock); verification_seconds=verification_seconds+check_seconds;
                returned_rms=norm(candidate.phi-current.phi)/sqrt(N^3);
                fallback(branch)=lb_kick_can_increase(false,candidate.gate.verified,candidate.index, ...
                    source,returned_rms,step_kick,candidate.energy,current.energy);
                row=lb_trial_record('lbsd4_down',label,target,sign,candidate,toc(attempt_clock),check_seconds);
                row.initial_rms=initial_rms; row.kick_rms=kick_rms; row.initial_drive=initial_drive;
                row.return_parent_rms=returned_rms; row.gate_status=string(candidate.gate.status);
                attempt_rows{end+1}=row; %#ok<SAGROW>
                lb_write_stage(ctx.run_dir,['lbsd_' label],candidate,row);
                if candidate.valid
                    current=candidate; accepted=true; break;
                end
            end
            if accepted, break; end
            if ~all(fallback)
                fprintf('NO_KICK_INCREASE: unresolved or non-parent branch at amplitude %.4g\n',step_kick);
                break;
            end
        end
        if ~accepted, current.valid=false; break; end
        stages{end+1}=current.phi; stage_indices(end+1)=current.index; %#ok<SAGROW>
    end
    sequential_complete=current.valid && isequal(stage_indices,[4 3 2 1]);
    current.valid=sequential_complete;
    search_seconds=toc(method_clock);

    %% Audit the saved raw negative mode and integrate both fresh GF branches.
    endpoints=struct([]); negative_verified=false;
    if sequential_complete
        v1=current.V_full(:,1); v1=v1-mean(v1); v1=v1/(norm(v1)/sqrt(N^3));
        check_clock=tic; Hv=lb_hessian_action(current.phi,v1);
        direction_rayleigh=real(v1'*Hv)/N^3;
        direction_residual=norm(Hv-direction_rayleigh*v1)/sqrt(N^3);
        negative_verified=direction_rayleigh<-1e-7 && direction_residual<=1e-8;
        verification_seconds=verification_seconds+toc(check_clock);
    elseif ctx.quick
        v1=current.out.V(:,1); v1=v1-mean(v1); v1=v1/(norm(v1)/sqrt(N^3));
    end
    current.negative_mode_verified=negative_verified;
    gf_seconds=0;
    if (sequential_complete && negative_verified) || ctx.quick
        gf=struct('dt',1,'maxiter',200000,'tol',1e-9, ...
            'energy_every',100,'flush_every',200,'checkpoint_every',2000,'nesterov',false);
        if ctx.quick, gf.maxiter=3; end
        for branch=1:2
            branch_signs=[1 -1]; sign=branch_signs(branch);
            phi_start=current.phi+sign*1e-5*v1;
            branch_clock=tic;
            try
                [phi_end,flow]=lb_hex_flow_with_checkpoints(phi_start,gf,ctx.run_dir, ...
                    sprintf('lbsd_gf_sign_%+d',sign),current.phi,v1,sign);
            catch exception
                phi_end=phi_start;
                flow=struct('converged',false,'status','failed','error',exception.message, ...
                    'elapsed_seconds',toc(branch_clock),'options',gf, ...
                    'parent_phi_sha256',lb_numeric_sha(current.phi),'branch_sign',sign);
                fprintf('GF sign=%+d failed: %s\n',sign,exception.message);
            end
            gf_seconds=gf_seconds+flow.elapsed_seconds; check_clock=tic;
            endpoint=lb_branch_info(phi_end,seed.phi_bcc,current.phi,sign,hex_reference);
            endpoint.converged=endpoint.converged && flow.converged;
            verification_seconds=verification_seconds+toc(check_clock);
            endpoint.flow=flow;
            if isempty(endpoints), endpoints=endpoint; else, endpoints(branch)=endpoint; end
        end
    end
    routes=lb_make_route('bcc2hex','lbsd4_down',seed,current,stages,endpoints,ctx.quick,stage_indices);
    lbsd_rows=vertcat(attempt_rows{:});
    timing_rows{1}=struct('method',"lbsd4_down",'version',"hex1490068_exact_V_sequential", ...
        'iterations',sum([lbsd_rows.iterations]),'search_seconds',search_seconds, ...
        'kernel_seconds',sum([lbsd_rows.kernel_seconds]),'verification_seconds',verification_seconds, ...
        'gf_seconds',gf_seconds,'elapsed_before_plot_seconds',toc(method_clock),'complete',routes.complete);

    %% FFT-cached SR-HiSD 1490067: exact input, current tangent, one Ritz step.
    if ctx.quick
        seed=lb_demo_initial('bcc2hex','srhisd1',kick);
        V_seed=seed.V_bcc_tr;
    else
        stored=load(fullfile(ctx.data_dir,'input','bcc2hex_srhisd.mat'),'input'); input=stored.input;
        assert(input.N==128 && abs(input.L-24*sqrt(2)*pi)<1e-12);
        assert(input.par.tau==-.3 && input.par.gamma==.8);
        seed=struct('N',input.N,'L',input.L,'par',input.par,'kick',.01, ...
            'phi_bcc',input.phi_parent(:),'phi_init',input.phi_init(:), ...
            'bb',struct('dt',.24,'betau',.2/.24,'betat',2/.24, ...
            'maxiter',80000,'epsf',1e-7,'outputp',100));
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

    %% Recompute the raw-H negative direction before starting fresh GF.
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
    row.initial_rms=NaN; row.kick_rms=NaN; row.initial_drive=NaN;
    row.return_parent_rms=NaN; row.gate_status="raw_negative_check";
    attempt_rows{end+1}=row;
    lb_write_stage(ctx.run_dir,'srhisd1_final',current,row);
    save(fullfile(ctx.run_dir,'srhisd1_negative_mode.mat'),'spectrum','verification_seconds','-v7.3');

    endpoints=struct([]); gf_seconds=0;
    if current.valid || ctx.quick
        if negative_verified, v1=spectrum.V(:,1); else, v1=solver.V(:,1); end
        innp=@(X,Y) (X'*Y)/N^3; v1=lb_normalize_rms(v1,innp);
        gf=struct('dt',1,'maxiter',200000,'tol',1e-9,'outputEvery',2000, ...
            'store_history',false,'nesterov',false,'chunk_iter',2000);
        if ctx.quick, gf.maxiter=3; end
        drawing=struct('draw_fig',false,'verbose',true);
        for branch=1:2
            branch_signs=[1 -1]; sign=branch_signs(branch);
            [phi_end,flow]=lb_flow_with_checkpoints(phi+sign*1e-5*v1,gf,drawing, ...
                ctx.run_dir,sprintf('srhisd1_gf_sign_%+d',sign),phi,v1,sign);
            gf_seconds=gf_seconds+flow.elapsed_seconds; check_clock=tic;
            endpoint=lb_branch_info(phi_end,seed.phi_bcc,phi,sign,hex_reference);
            verification_seconds=verification_seconds+toc(check_clock);
            endpoint.flow=flow;
            if isempty(endpoints), endpoints=endpoint; else, endpoints(branch)=endpoint; end
        end
    end
    routes(2)=lb_make_route('bcc2hex','srhisd1',seed,current,{phi},endpoints,ctx.quick,1);
    timing_rows{2}=struct('method',"srhisd1",'version',"sr1490067_fft_cache", ...
        'iterations',solver.it,'search_seconds',solver.search_seconds,'kernel_seconds',solver.kernel_seconds, ...
        'verification_seconds',verification_seconds,'gf_seconds',gf_seconds, ...
        'elapsed_before_plot_seconds',toc(method_clock),'complete',routes(2).complete);
end

%% Save actual results, timings and automatically draw this case.
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
review_result=review_finish(ctx,struct('status',status,'completed',completed, ...
    'route_count',numel(routes),'numerical_search_performed',strcmp(ctx.mode,'compute'), ...
    'parameter_version','hex1490068_20261120','lbsd_stage_sequence',[4 3 2 1], ...
    'sr_implementation','sr1490067_fft_cache','figure_files',{figure_files}));
