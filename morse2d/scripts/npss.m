%% Retained-history NPSS comparison on the same Table 1 initial states
% Two thresholds; Stage I refreshes the spectral nullspace, Stage II uses H.
if ~exist('review_options','var'), review_options=struct(); end
task_dir=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(task_dir),'common'));
ctx=review_setup(mfilename('fullpath'),'npss',review_options);

%% Model and two-stage NPSS settings
c=struct('a',6,'r0',1,'dt',0.02,'maxiter',50000,'tol',1e-8, ...
    'res_tol',1e-8,'grad_tol',1e-5,'eig_neg_tol',1e-6,'cluster_cutoff',1.6, ...
    'bb_min_factor',0.001,'bb_max_factor',2,'null_tol',1e-6, ...
    'neg_tol',1e-6,'rayleigh_tol',0.1,'angle_sin_tol',0.25);
thresholds=[0.10 0.01];
Ns=5:10; n_trials=500;
if ctx.quick, Ns=5; n_trials=2; c.maxiter=50; end
run_info=struct('subtask','npss','mode',ctx.mode,'quick',ctx.quick,'Ns',Ns, ...
    'trials_per_N',n_trials,'thresholds',thresholds,'actual_solver_runs',0);
for N=Ns
    initial=load(fullfile(ctx.data_dir,'initial_states',sprintf('n%d.mat',N)));
    assert(isequal(size(initial.q0_list),[2*N 500]));
end
saved_file=fullfile(ctx.data_dir,'reference','npss','result.mat');
if ~strcmp(ctx.mode,'compute'), assert(isfile(saved_file)); end
if strcmp(ctx.mode,'check')
    review_result=review_finish(ctx,run_info); return;
end

%% Load the current full calculation when no new searches are requested
if strcmp(ctx.mode,'reference')
    saved=load(saved_file,'Trials','Rates','Endpoints','c','run_info');
    Trials=saved.Trials; Rates=saved.Rates; Endpoints=saved.Endpoints; c=saved.c;
    run_info.source='saved_full_calculation';
    run_info.source_file=saved_file;
    run_info.source_solver_runs=saved.run_info.actual_solver_runs;
    save(fullfile(ctx.run_dir,'result.mat'),'Trials','Rates','Endpoints','c','run_info','-v7');
    writetable(Trials,fullfile(ctx.run_dir,'trials.csv'));
    writetable(Rates,fullfile(ctx.run_dir,'rates.csv'));
    disp(Rates);
    review_result=review_finish(ctx,run_info); return;
end

%% Run each threshold on all shared initial states; certify the actual endpoints
rows={}; Endpoints=cell(2,numel(Ns));
for threshold_id=1:2
    c.rayleigh_tol=thresholds(threshold_id);
    method=string(sprintf('NPSS, epsilon_lambda=%.2f',c.rayleigh_tol));
    for iN=1:numel(Ns)
        N=Ns(iN);
        initial=load(fullfile(ctx.data_dir,'initial_states',sprintf('n%d.mat',N)));
        endpoints=zeros(2*N,n_trials);
        for trial=1:n_trials
            started=tic;
            [q,~,~,iterations,info]=morse_npss(initial.q0_list(:,trial),initial.qmin,c);
            elapsed=toc(started);
            cert=morse_certificate(q,initial.qmin,c);
            assert(info.bb_history_resets==0,'Retained-history NPSS must not reset BB secants.');
            rows{end+1}=struct('method',method,'epsilon_lambda',c.rayleigh_tol,'N',N, ...
                'trial',trial,'success',cert.success,'mother',cert.mother, ...
                'converged',cert.converged,'rms',cert.rms,'gn',cert.gn,'index',cert.index, ...
                'iterations',iterations,'stage',info.stage,'null_updates',info.null_updates, ...
                'bb_history_resets',info.bb_history_resets,'seconds',elapsed); %#ok<SAGROW>
            endpoints(:,trial)=q;
            if mod(trial,100)==0 || ctx.quick
                fprintf('%s N=%d trial=%d/%d RMS=%.3e index=%d hit=%d mother=%d\n', ...
                    method,N,trial,n_trials,cert.rms,cert.index,cert.success,cert.mother);
            end
        end
        Endpoints{threshold_id,iN}=endpoints;
        save(fullfile(ctx.run_dir,'checkpoint.mat'),'rows','Endpoints','c', ...
            'run_info','threshold_id','iN','-v7');
    end
end

%% Count hits, mother returns and unresolved trials separately
Trials=struct2table(vertcat(rows{:})); summary_rows={};
for threshold_id=1:2
    for N=Ns
        group=Trials(Trials.epsilon_lambda==thresholds(threshold_id) & Trials.N==N,:);
        summary_rows{end+1}=struct('method',group.method(1),'N',N,'trials',height(group), ...
            'successes',sum(group.success),'success_pct',100*mean(group.success), ...
            'mother_returns',sum(group.mother),'unconverged',sum(~group.converged)); %#ok<SAGROW>
    end
end
Rates=struct2table(vertcat(summary_rows{:}));
run_info.actual_solver_runs=strcmp(ctx.mode,'compute')*height(Trials);
save(fullfile(ctx.run_dir,'result.mat'),'Trials','Rates','Endpoints','c','run_info','-v7');
writetable(Trials,fullfile(ctx.run_dir,'trials.csv'));
writetable(Rates,fullfile(ctx.run_dir,'rates.csv'));
disp(Rates);
review_result=review_finish(ctx,run_info);
