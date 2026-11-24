function [phi,flow]=lb_flow_with_checkpoints(phi_start,options,drawing,output_dir,tag,parent_phi,parent_direction,branch_sign)
%LB_FLOW_WITH_CHECKPOINTS Run fresh GF chunks and atomically save each chunk.
% Numerical updates are performed only by lb_gradient_flow. Its last chunk
% may stop early without reporting the exact update count, so bounds are
% retained instead of inventing a last-chunk iteration count. This function
% does not resume an earlier checkpoint.
global LB
assert(all(isfield(options,{'dt','maxiter','tol'})), ...
    'LB:FlowOptions','GF options must explicitly contain dt, maxiter and tol.');
assert(isscalar(options.dt) && options.dt>0 && isfinite(options.dt) && ...
    isscalar(options.maxiter) && options.maxiter>=1 && options.maxiter==fix(options.maxiter) && ...
    isscalar(options.tol) && options.tol>0 && isfinite(options.tol), ...
    'LB:FlowOptions','Invalid GF step, iteration budget or tolerance.');
assert(~isfield(options,'nesterov') || ~options.nesterov, ...
    'LB:FlowMomentum','Chunked GF uses the supplied non-momentum protocol.');
assert(isscalar(branch_sign) && ismember(branch_sign,[-1 1]), ...
    'LB:FlowBranch','The branch sign must be +1 or -1.');
tag=char(string(tag));
assert(~isempty(regexp(tag,'^[A-Za-z0-9_+-]+$','once')), ...
    'LB:FlowTag','Use a simple file tag with letters, numbers, underscores and signs.');
assert(isfolder(output_dir),'LB:FlowOutput','The task output directory must exist.');
phi=phi_start(:); parent_phi=parent_phi(:); parent_direction=parent_direction(:);
assert(numel(phi)==LB.N^3 && isequal(size(phi),size(parent_phi),size(parent_direction)) && ...
    all(isfinite(phi)) && all(isfinite(parent_phi)) && all(isfinite(parent_direction)), ...
    'LB:FlowState','The start, parent and direction must be finite grid vectors.');
chunk_iter=2000;
if isfield(options,'chunk_iter'), chunk_iter=options.chunk_iter; end
assert(isscalar(chunk_iter) && chunk_iter>=1 && chunk_iter==fix(chunk_iter), ...
    'LB:FlowChunk','The checkpoint chunk size must be a positive integer.');
destination=fullfile(output_dir,[tag '_checkpoint.mat']);
assert(~isfile(destination),'LB:FlowExistingCheckpoint', ...
    'Fresh GF refuses an existing checkpoint; choose a new run: %s',destination);
temporary=fullfile(output_dir,[tag '_checkpoint.tmp.mat']);
parent_hash=lb_numeric_sha(parent_phi);
direction_hash=lb_numeric_sha(parent_direction);
initial_phi_hash=lb_numeric_sha(phi);
gf_options=options; gf_options.chunk_iter=chunk_iter;
N=LB.N; L=LB.L; par=LB.par;
lower=0; upper=0;
history=struct('updates_lower_bound',{},'updates_upper_bound',{}, ...
    'force_rms',{},'energy',{},'elapsed_seconds',{});
clock=tic;
while upper<options.maxiter
    force=lb_negative_gradient(phi); residual=sqrt(mean(force.^2));
    assert(all(isfinite(phi)) && isfinite(residual),'LB:FlowNonfinite','Nonfinite GF state.');
    if residual<options.tol, break; end
    capacity=min(chunk_iter,options.maxiter-upper);
    chunk_options=options; chunk_options.maxiter=capacity;
    chunk_options.tag=sprintf('%s-after%d',tag,upper);
    [phi,~]=lb_gradient_flow(phi,chunk_options,drawing);
    force=lb_negative_gradient(phi); residual=sqrt(mean(force.^2));
    assert(all(isfinite(phi)) && isfinite(residual),'LB:FlowNonfinite','Nonfinite evolved GF state.');
    if residual<options.tol
        upper=upper+capacity;
    else
        lower=lower+capacity; upper=lower;
    end
    elapsed_seconds=toc(clock);
    history(end+1)=struct('updates_lower_bound',lower,'updates_upper_bound',upper, ...
        'force_rms',residual,'energy',lb_energy(phi),'elapsed_seconds',elapsed_seconds); %#ok<AGROW>
    write_checkpoint();
    fprintf('GF CHECKPOINT %s sign%+d bounds=[%d,%d] RMS=%.6e E=%.12g seconds=%.1f\n', ...
        tag,branch_sign,lower,upper,residual,history(end).energy,elapsed_seconds);
end
force=lb_negative_gradient(phi); residual=sqrt(mean(force.^2));
elapsed_seconds=toc(clock);
converged=all(isfinite(phi)) && isfinite(residual) && residual<options.tol;
write_checkpoint();
flow=struct('converged',converged,'elapsed_seconds',elapsed_seconds, ...
    'updates_lower_bound',lower,'updates_upper_bound',upper, ...
    'update_count_exact',lower==upper,'force_rms',residual,'energy',lb_energy(phi), ...
    'history',history,'parent_phi_sha256',parent_hash,'direction_sha256',direction_hash, ...
    'initial_phi_sha256',initial_phi_hash,'branch_sign',branch_sign,'options',gf_options, ...
    'checkpoint_file',destination,'resumed',false);

    function write_checkpoint()
        % Same-directory rename keeps a complete previous checkpoint until
        % the new MAT has been completely written.
        checkpoint=struct('phi',phi,'lower',lower,'upper',upper,'history',history, ...
            'elapsed_seconds',elapsed_seconds,'parent_hash',parent_hash, ...
            'direction_hash',direction_hash,'initial_phi_hash',initial_phi_hash, ...
            'branch_sign',branch_sign,'gf_options',gf_options,'N',N,'L',L,'par',par);
        save(temporary,'-struct','checkpoint','-v7.3');
        movefile(temporary,destination,'f');
    end
end
