%% Five original Table 1 methods: HiSD orders 1--4 and one SR-HiSD direction
% Run this script directly, or select morse2d/original_methods in main.m.
% The mathematical experiment is explicit below; setup only handles paths/logs.
if ~exist('review_options','var'), review_options=struct(); end
task_dir=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(task_dir),'common'));
ctx=review_setup(mfilename('fullpath'),'original_methods',review_options);

%% Model, perturbations and stopping tests
c=struct('a',6,'r0',1,'dt',0.02,'maxiter',50000,'tol',1e-8, ...
    'res_tol',1e-8,'grad_tol',1e-5,'eig_neg_tol',1e-6,'cluster_cutoff',1.6, ...
    'bb_min_factor',0.001,'bb_max_factor',2);
Ns=5:10; n_trials=500; random_seed=20260714; perturb_radius=0.1;
if ctx.quick, Ns=5; n_trials=2; c.maxiter=50; end
methods=["1-SR-HiSD","1-HiSD","2-HiSD","3-HiSD","4-HiSD"];
direction_counts=[1 1 2 3 4];
run_info=struct('subtask','original_methods','mode',ctx.mode,'quick',ctx.quick, ...
    'Ns',Ns,'trials_per_N',n_trials,'methods',methods,'actual_solver_runs',0);

%% Input check: every method uses the same frozen mother and initial states
for N=Ns
    inputs=load(fullfile(ctx.data_dir,'initial_states',sprintf('n%d.mat',N)));
    assert(isequal(size(inputs.q0_list),[2*N 500]) && numel(inputs.qmin)==2*N);
    assert(isfile(fullfile(ctx.data_dir,'minima',sprintf('n%d.mat',N))));
end
saved_file=fullfile(ctx.data_dir,'reference','original_methods','result.mat');
if ~strcmp(ctx.mode,'compute'), assert(isfile(saved_file)); end
if strcmp(ctx.mode,'check')
    review_result=review_finish(ctx,run_info); return;
end

%% Load the current full calculation when only saved results are requested
if strcmp(ctx.mode,'reference')
    saved=load(saved_file,'Trials','Rates','InitialStates','Endpoints','c','run_info');
    Trials=saved.Trials; Rates=saved.Rates; InitialStates=saved.InitialStates;
    Endpoints=saved.Endpoints; c=saved.c;
    run_info.source='saved_full_calculation';
    run_info.source_file=saved_file;
    run_info.source_solver_runs=saved.run_info.actual_solver_runs;
    save(fullfile(ctx.run_dir,'result.mat'),'Trials','Rates','InitialStates','Endpoints','c','run_info','-v7');
    writetable(Trials,fullfile(ctx.run_dir,'trials.csv'));
    writetable(Rates,fullfile(ctx.run_dir,'rates.csv'));
    disp(Rates);
    review_result=review_finish(ctx,run_info); return;
end

%% Relax one mother minimum, sample 500 quotient perturbations, run five methods
rows={}; InitialStates=cell(numel(Ns),1); Endpoints=cell(numel(Ns),1);
for iN=1:numel(Ns)
    N=Ns(iN);
    shared=load(fullfile(ctx.data_dir,'initial_states',sprintf('n%d.mat',N)));
    minimum=load(fullfile(ctx.data_dir,'minima',sprintf('n%d.mat',N)),'global_min_q');
    [qmin,~,~]=morse_gradient_descent(minimum.global_min_q(:),5e-4,12000,c.a,c.r0,c.tol);
    [~,~,Qmin]=morse_rigid_basis(qmin);
    rng(random_seed,'twister');
    initial_states=zeros(2*N,n_trials); endpoints=zeros(2*N,n_trials,5);
    for trial=1:n_trials
        % An isotropic fixed-radius perturbation after removing rigid modes.
        z=randn(size(Qmin,2),1); z=z/norm(z);
        q0=qmin+perturb_radius*(Qmin*z);
        positions=reshape(q0,2,N)'; positions=positions-mean(positions,1);
        q0=reshape(positions',[],1);
        assert(isequal(q0,shared.q0_list(:,trial)), ...
            'Review:MorseInitialState','Generated initial state differs from the shared input.');
        initial_states(:,trial)=q0;
        for method_id=1:numel(methods)
            % HiSD starts in the full Hessian space; SR-HiSD in the quotient.
            [~,~,H0]=morse_energy_gradient_hessian(q0,c.a,c.r0);
            if method_id==1
                [~,~,Q0]=morse_rigid_basis(q0);
                [Z,D]=eig(Q0'*H0*Q0); [~,order]=sort(diag(D),'ascend');
                V0=Q0*Z(:,order); V0=V0(:,1);
            else
                [Z,D]=eig(H0); [~,order]=sort(diag(D),'ascend');
                V0=Z(:,order(1:direction_counts(method_id)));
            end
            fprintf('N=%d trial=%d/%d method=%s\n',N,trial,n_trials,methods(method_id));
            started=tic;
            if method_id==1
                [q,~,residual,iterations]=morse_srhisd(q0,V0,c.dt,c.maxiter,c.a,c.r0,c);
            else
                [q,~,residual,iterations]=morse_hisd(q0,V0,c.dt,c.maxiter,c.a,c.r0,c);
            end
            elapsed=toc(started);
            cert=morse_certificate(q,qmin,c);
            assert(abs(residual-cert.rms)<=1e-12,'Returned residual must describe the endpoint.');
            rows{end+1}=struct('method',methods(method_id),'N',N,'trial',trial, ...
                'success',cert.success,'rms',cert.rms,'gn',cert.gn,'index',cert.index, ...
                'connected',cert.components==1,'iterations',iterations,'seconds',elapsed); %#ok<SAGROW>
            endpoints(:,trial,method_id)=q;
            fprintf('  index=%d RMS=%.3e index-one hit=%d\n',cert.index,cert.rms,cert.success);
        end
    end
    InitialStates{iN}=initial_states; Endpoints{iN}=endpoints;
    save(fullfile(ctx.run_dir,'checkpoint.mat'),'rows','InitialStates','Endpoints', ...
        'c','run_info','iN','-v7');
end

%% Aggregate every trial, preserving failures in the denominator
Trials=struct2table(vertcat(rows{:})); summary_rows={};
for method_id=1:numel(methods)
    for N=Ns
        group=Trials(Trials.method==methods(method_id) & Trials.N==N,:);
        summary_rows{end+1}=struct('method',methods(method_id),'N',N,'trials',height(group), ...
            'successes',sum(group.success),'success_pct',100*mean(group.success)); %#ok<SAGROW>
    end
end
Rates=struct2table(vertcat(summary_rows{:}));
run_info.actual_solver_runs=strcmp(ctx.mode,'compute')*height(Trials);
save(fullfile(ctx.run_dir,'result.mat'),'Trials','Rates','InitialStates','Endpoints','c','run_info','-v7');
writetable(Trials,fullfile(ctx.run_dir,'trials.csv'));
writetable(Rates,fullfile(ctx.run_dir,'rates.csv'));
disp(Rates);
review_result=review_finish(ctx,run_info);
