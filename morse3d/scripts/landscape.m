%% Build the N=8 Morse landscape: comparison seeds -> upward search -> descent.
% The queue, solver calls, acceptance tests and merging rules are visible here.
if ~exist('review_options','var'); review_options=struct(); end
task_dir=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(task_dir),'common'));
ctx=review_setup(mfilename('fullpath'),'landscape',review_options);
paper_validation=[];

%% Published settings (Quick changes only the execution budget).
N=8; a=6; r0=1; tolerance=1e-8; dt=5e-4;
maxiter=60000; max_depth=3; max_queue=28; max_index=3;
seed_perturbation=5e-7; signs=[1 -1]; unique_tol=5e-2;
up_perturbations=[5e-8 1e-7 5e-7 1e-6 5e-6 1e-5 5e-5];
descent_steps=30000; descent_dt=[5e-4 2e-4 1e-4];
descent_perturbations=[1e-4 5e-4 1e-3 2e-3];
if ctx.quick
    maxiter=10; max_depth=1; max_queue=1; max_index=1;
    signs=1; up_perturbations=5e-7; descent_steps=10;
end
minimum=load(fullfile(ctx.data_dir,'input','minimum.mat'),'global_min_q');
assert(numel(minimum.global_min_q)==3*N && all(isfinite(minimum.global_min_q(:))));
assert(isfile(fullfile(ctx.data_dir,'paper','figure1_states.csv')));
if strcmp(ctx.mode,'check')
    review_result=review_finish(ctx,struct('status','PASS','paper_states',17,'paper_edges',106));
    return;
end
if strcmp(ctx.mode,'reference')
    saved=load(fullfile(ctx.data_dir,'reference','result.mat'),'Nodes','Edges','CompareRows','qmin');
    Nodes=saved.Nodes; Edges=saved.Edges; CompareRows=saved.CompareRows;
    qmin=saved.qmin; numerical_search_performed=false;
else
    %% Reuse this run's comparison seeds, or generate them when run on its own.
    seed_file=review_find_result(ctx.task_dir,'compare_methods',ctx.options);
    if ~isempty(seed_file)
        seeds=load(seed_file,'Nodes','Edges','CompareRows','qmin','settings','numerical_search_performed');
        assert(isfield(seeds.settings,'perturbation') && seeds.numerical_search_performed && seeds.settings.quick==ctx.quick, ...
            'Morse3D:SeedMode','Comparison seeds must be a matching fresh computation.');
        Nodes=seeds.Nodes; Edges=seeds.Edges; CompareRows=seeds.CompareRows; qmin=seeds.qmin;
        fprintf('Continuing the comparison graph from %s\n',seed_file);
    else
        qmin=minimum.global_min_q;
        for attempt=0:3
            [qmin,min_res,min_iter]=morse3d_gradient_descent(qmin,5e-4/(2^attempt), ...
                12000*(1+(attempt>0)),a,r0,tolerance);
            [Emin,gmin,Hmin]=morse3d_energy(qmin,a,r0);
            [~,~,Qmin]=morse3d_rigid_basis(qmin);
            spectrum=sort(eig((Qmin'*Hmin*Qmin+(Qmin'*Hmin*Qmin)')/2));
            root_ok=all(isfinite([Emin;gmin(:);min_res])) && ...
                sum(spectrum < -1e-6)==0 && norm(gmin)<1e-5 && min_res<5e-8;
            if root_ok; break; end
        end
        assert(root_ok,'Morse3D:Minimum','Root minimum failed strict convergence.');
        root_node=morse3d_make_node(qmin,Emin,norm(gmin),0,spectrum,'minimizer','minimum', ...
            0,0,0,5e-4,0,1,min_iter,min_res,true,true,false,false);
        [Nodes,root_id]=morse3d_merge(struct([]),root_node,unique_tol);
        Edges=struct([]); CompareRows=struct([]);
        Vminimum=morse3d_initial_frame(qmin,a,r0,true);
        for sign_value=signs
            q0=qmin+sign_value*seed_perturbation*Vminimum(:,1);
            for method={'SR-HiSD','HiSD'}
                algorithm=method{1}; started=tic;
                if strcmp(algorithm,'SR-HiSD')
                    V0=morse3d_initial_frame(q0,a,r0,true);
                    [q,V,res,iterations]=morse3d_srhisd(q0,V0(:,1),dt,maxiter,a,r0,tolerance);
                else
                    V0=morse3d_initial_frame(q0,a,r0,false);
                    [q,V,res,iterations]=morse3d_hisd(q0,V0(:,1:4),dt,maxiter,a,r0,tolerance);
                end
                [q,res,iterations]=morse3d_continue(algorithm,q,dt,maxiter,iterations,res,1,a,r0,tolerance);
                row=morse3d_record(q,q0,qmin,a,r0,res,iterations,algorithm,dt, ...
                    seed_perturbation,sign_value,1,1);
                row.time_sec=toc(started);
                [Nodes,id,is_new]=morse3d_merge(Nodes,row.node,unique_tol);
                row.source_id=root_id; row.target_node_id=id; row.target_is_new=is_new;
                if isempty(CompareRows); CompareRows=row; else; CompareRows(end+1)=row; end
                Edges=morse3d_add_edge(Edges,Nodes,root_id,id,row,1,1);
            end
        end
    end

    %% Upward queue: search index k+1 from every accepted index-k state.
    root_id=1; queue=root_id; queued=false(numel(Nodes),1); queued(root_id)=true;
    % Include accepted SR comparison targets exactly as in the original search.
    for j=1:numel(CompareRows)
        row=CompareRows(j);
        if strcmp(row.algo,'SR-HiSD') && row.target_ok && ~queued(row.target_node_id)
            queue(end+1)=row.target_node_id; queued(row.target_node_id)=true;
        end
    end
    pointer=1; expanded=0;
    while pointer<=numel(queue) && expanded<max_queue
        parent_id=queue(pointer); pointer=pointer+1; expanded=expanded+1;
        parent=Nodes(parent_id);
        if parent.depth>=max_depth; continue; end
        Vparent=morse3d_initial_frame(parent.q,a,r0,true);
        k=min(max(parent.morse_index+1,1),size(Vparent,2));
        if k>max_index; continue; end
        for amplitude=up_perturbations
            for sign_value=signs
                q0=parent.q+sign_value*amplitude*Vparent(:,k);
                V0=morse3d_initial_frame(q0,a,r0,true); started=tic;
                [q,V,res,iterations]=morse3d_srhisd(q0,V0(:,1:k),dt,maxiter,a,r0,tolerance);
                [q,res,iterations]=morse3d_continue('SR-HiSD',q,dt,maxiter,iterations,res,k,a,r0,tolerance);
                row=morse3d_record(q,q0,qmin,a,r0,res,iterations,'SR-HiSD',dt, ...
                    amplitude,sign_value,parent.depth+1,k);
                row.time_sec=toc(started); row.k_target=k; row.parent_id=parent_id;
                [Nodes,id,is_new]=morse3d_merge(Nodes,row.node,unique_tol);
                accepted=row.stable_ok && row.morse_index==k && row.gn<1e-5 && row.res<5e-8;
                row.target_ok=accepted; row.target_node_id=id; row.target_is_new=is_new;
                Nodes(id).target_ok=Nodes(id).target_ok || accepted;
                CompareRows(end+1)=row;
                Edges=morse3d_add_edge(Edges,Nodes,parent_id,id,row,is_new,parent.depth+1);
                if accepted && is_new
                    if id>numel(queued) || ~queued(id)
                        queue(end+1)=id; queued(id)=true;
                    end
                end
                fprintf('parent=%d -> index %d: E=%.12g, accepted=%d\n',parent_id,k,row.E,accepted);
            end
        end
    end

    %% Bilateral gradient flow from every accepted index-one saddle.
    saddle_ids=find([Nodes.morse_index]==1 & [Nodes.target_ok]);
    for saddle_id=saddle_ids
        saddle=Nodes(saddle_id); Vunstable=morse3d_initial_frame(saddle.q,a,r0,true);
        for sign_value=[-1 1]
            found=false; best_energy=inf; best_row=struct();
            for amplitude=descent_perturbations
                q0=saddle.q+sign_value*amplitude*Vunstable(:,1);
                for flow_dt=descent_dt
                    started=tic;
                    [q,res,iterations]=morse3d_gradient_descent(q0,flow_dt,descent_steps,a,r0,tolerance);
                    if iterations>=descent_steps && isfinite(res) && res<5e-3
                        [q,res,extra]=morse3d_gradient_descent(q,.5*flow_dt,2*descent_steps,a,r0,tolerance);
                        iterations=iterations+extra;
                    end
                    row=morse3d_record(q,q0,qmin,a,r0,res,iterations,'GD',flow_dt, ...
                        amplitude,sign_value,saddle.depth+1,0);
                    row.time_sec=toc(started);
                    if row.target_ok && row.E<best_energy
                        best_row=row; best_energy=row.E; found=true;
                    end
                end
            end
            if ~found; continue; end
            row=best_row; row.parent_id=saddle_id; row.k_target=0;
            [Nodes,id,is_new]=morse3d_merge(Nodes,row.node,unique_tol,true,r0);
            row.target_node_id=id; row.target_is_new=is_new;
            CompareRows(end+1)=row;
            Edges=morse3d_add_edge(Edges,Nodes,saddle_id,id,row,is_new,saddle.depth+1);
        end
    end
    numerical_search_performed=true;
end
if ~ctx.quick
    % Check the article's displayed states and the calculated graph itself.
    paper_validation=morse3d_check_paper_result(Nodes,Edges,CompareRows);
end

%% Store the actual graph; plotting is an independent script.
settings=struct('N',N,'a',a,'r0',r0,'dt',dt,'maxiter',maxiter,'max_depth',max_depth, ...
    'max_queue',max_queue,'max_index',max_index,'up_perturbations',up_perturbations, ...
    'descent_dt',descent_dt,'descent_perturbations',descent_perturbations,'quick',ctx.quick);
result_file=morse3d_save_result(ctx.run_dir,struct('Nodes',Nodes,'Edges',Edges, ...
    'CompareRows',CompareRows,'qmin',qmin,'a',a,'r0',r0,'N',N, ...
    'settings',settings,'numerical_search_performed',numerical_search_performed));
if ~isempty(paper_validation)
    save(result_file,'paper_validation','-append');
    assert(paper_validation.passed,'Morse3D:PaperMismatch', ...
        'Graph does not satisfy the article checks; calculated states are retained.');
end
review_result=review_finish(ctx,struct('status','PASS','result_file',result_file, ...
    'nodes',numel(Nodes),'edges',numel(Edges),'trials',numel(CompareRows), ...
    'numerical_search_performed',numerical_search_performed));
