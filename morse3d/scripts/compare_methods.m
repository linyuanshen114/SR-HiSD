%% Compare full-space HiSD and SR-HiSD from identical N=8 perturbations.
% Open this script and press Run, or select morse3d/compare_methods in main.m.
if ~exist('review_options','var'); review_options=struct(); end
task_dir=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(task_dir),'common'));
ctx=review_setup(mfilename('fullpath'),'compare_methods',review_options);

%% Model and numerical settings.
N=8; a=6; r0=1; tolerance=1e-8;
dt=5e-4; perturbation=5e-7; signs=[1 -1]; maxiter=60000;
if ctx.quick; signs=1; maxiter=10; end
minimum=load(fullfile(ctx.data_dir,'input','minimum.mat'),'global_min_q');
assert(numel(minimum.global_min_q)==3*N && all(isfinite(minimum.global_min_q(:))));
if strcmp(ctx.mode,'check')
    review_result=review_finish(ctx,struct('status','PASS','grid_dof',3*N));
    return;
end

if strcmp(ctx.mode,'reference')
    % Replay the current validated four-record comparison result.
    saved=load(fullfile(ctx.data_dir,'reference','compare_methods','result.mat'), ...
        'Nodes','Edges','CompareRows','qmin');
    Nodes=saved.Nodes; Edges=saved.Edges; CompareRows=saved.CompareRows; qmin=saved.qmin;
    numerical_search_performed=false;
else
    %% Relax the supplied minimum, then remove translations and rotations.
    qmin=minimum.global_min_q;
    for attempt=0:3
        relax_dt=5e-4/(2^attempt); relax_steps=12000*(1+(attempt>0));
        [qmin,min_res,min_iter]=morse3d_gradient_descent(qmin,relax_dt,relax_steps,a,r0,tolerance);
        [Emin,gmin,Hmin]=morse3d_energy(qmin,a,r0);
        [~,~,Qmin]=morse3d_rigid_basis(qmin);
        spectrum=sort(eig((Qmin'*Hmin*Qmin+(Qmin'*Hmin*Qmin)')/2));
        root_ok=all(isfinite([Emin;gmin(:);min_res])) && ...
            sum(spectrum < -1e-6)==0 && norm(gmin)<1e-5 && min_res<5e-8;
        if root_ok; break; end
    end
    assert(root_ok,'Morse3D:Minimum','The supplied minimum did not converge.');
    root_node=morse3d_make_node(qmin,Emin,norm(gmin),0,spectrum,'minimizer','minimum', ...
        0,0,0,5e-4,0,1,min_iter,min_res,true,true,false,false);
    [Nodes,root_id]=morse3d_merge(struct([]),root_node,5e-2);
    Edges=struct([]); CompareRows=struct([]);
    Vminimum=morse3d_initial_frame(qmin,a,r0,true);

    %% Both methods see exactly the same signed position perturbation.
    for sign_value=signs
        q0=qmin+sign_value*perturbation*Vminimum(:,1);
        for method={'SR-HiSD','HiSD'}
            algorithm=method{1}; started=tic;
            if strcmp(algorithm,'SR-HiSD')
                V0=morse3d_initial_frame(q0,a,r0,true);
                [q,V,res,iterations]=morse3d_srhisd(q0,V0(:,1),dt,maxiter,a,r0,tolerance);
            else
                V0=morse3d_initial_frame(q0,a,r0,false);
                [q,V,res,iterations]=morse3d_hisd(q0,V0(:,1:4),dt,maxiter,a,r0,tolerance);
            end
            % Original near-convergence continuation; certification uses new q.
            [q,res,iterations]=morse3d_continue(algorithm,q,dt,maxiter,iterations,res,1,a,r0,tolerance);
            row=morse3d_record(q,q0,qmin,a,r0,res,iterations,algorithm,dt, ...
                perturbation,sign_value,1,1);
            row.time_sec=toc(started);
            [Nodes,id,is_new]=morse3d_merge(Nodes,row.node,5e-2);
            row.source_id=root_id; row.target_node_id=id; row.target_is_new=is_new;
            if isempty(CompareRows); CompareRows=row; else; CompareRows(end+1)=row; end
            Edges=morse3d_add_edge(Edges,Nodes,root_id,id,row,1,1);
            fprintf('%s sign=%+d: E=%.12g, quotient index=%d, target=%d\n', ...
                algorithm,sign_value,row.E,row.morse_index,row.target_ok);
        end
    end
    numerical_search_performed=true;
end

%% Save numerical graph seeds for the independently runnable landscape step.
settings=struct('N',N,'a',a,'r0',r0,'dt',dt,'perturbation',perturbation, ...
    'signs',signs,'maxiter',maxiter,'tolerance',tolerance,'quick',ctx.quick);
result_file=morse3d_save_result(ctx.run_dir,struct('Nodes',Nodes,'Edges',Edges, ...
    'CompareRows',CompareRows,'qmin',qmin,'a',a,'r0',r0,'N',N, ...
    'settings',settings,'numerical_search_performed',numerical_search_performed));
review_result=review_finish(ctx,struct('status','PASS','result_file',result_file, ...
    'comparison_records',numel(CompareRows),'numerical_search_performed',numerical_search_performed));
