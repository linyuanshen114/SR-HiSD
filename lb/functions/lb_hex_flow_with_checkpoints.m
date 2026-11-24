function [phi,flow]=lb_hex_flow_with_checkpoints(phi_start,options,output_dir,tag,parent_phi,parent_direction,branch_sign)
%LB_HEX_FLOW_WITH_CHECKPOINTS Exact-update GF for the HEX sequential route.
% Keep the downloaded bilateral-flow update and stop checks, including its
% arithmetic order. The other LB cases retain their own GF implementation.
global LB
assert(all(isfield(options,{'dt','maxiter','tol'})), ...
    'LB:HexFlowOptions','Specify dt, maxiter, and tol explicitly.');
assert(isscalar(options.dt) && isfinite(options.dt) && options.dt>0 && ...
    isscalar(options.maxiter) && options.maxiter>=0 && options.maxiter==fix(options.maxiter) && ...
    isscalar(options.tol) && isfinite(options.tol) && options.tol>0, ...
    'LB:HexFlowOptions','Invalid time step, update budget, or force tolerance.');
options=with_defaults(options,struct('kick',1e-5,'energy_every',100, ...
    'flush_every',200,'checkpoint_every',2000,'nesterov',false));
assert(~options.nesterov,'LB:HexFlowMomentum','This GF protocol has no momentum.');
for key={'energy_every','flush_every','checkpoint_every'}
    value=options.(key{1});
    assert(isscalar(value) && isfinite(value) && value>=1 && value==fix(value), ...
        'LB:HexFlowOptions','Recording intervals must be positive integers.');
end
assert(isscalar(branch_sign) && ismember(branch_sign,[-1 1]), ...
    'LB:HexFlowBranch','The branch sign must be +1 or -1.');
tag=char(string(tag));
assert(~isempty(regexp(tag,'^[A-Za-z0-9_+-]+$','once')), ...
    'LB:HexFlowTag','Use a simple file tag.');
assert(isfolder(output_dir),'LB:HexFlowOutput','The output directory must exist.');
phi=phi_start(:); parent_phi=parent_phi(:); parent_direction=parent_direction(:);
m=numel(phi);
assert(m==LB.N^3 && isequal(size(phi),size(parent_phi),size(parent_direction)) && ...
    all(isfinite(phi)) && all(isfinite(parent_phi)) && all(isfinite(parent_direction)), ...
    'LB:HexFlowState','The start, parent, and direction must be finite grid vectors.');
history_file=fullfile(output_dir,[tag '_history.csv']);
checkpoint_file=fullfile(output_dir,[tag '_checkpoint.mat']);
endpoint_file=fullfile(output_dir,[tag '_endpoint.mat']);
initial_file=fullfile(output_dir,[tag '_initial.mat']);
assert(~any(cellfun(@isfile,{history_file,checkpoint_file,endpoint_file,initial_file})), ...
    'LB:HexFlowExisting','Fresh GF refuses existing branch output files.');
submitted_initial_hash=lb_numeric_sha(phi);
phi=phi-mean(phi);
parent_hash=lb_numeric_sha(parent_phi);
direction_hash=lb_numeric_sha(parent_direction);
initial_hash=lb_numeric_sha(phi);
dv=phi-parent_phi;
Fstart=lb_negative_gradient(phi); Fparent=lb_negative_gradient(parent_phi);
initial=struct('sign',branch_sign,'kick_RMS',norm(dv)/sqrt(m), ...
    'kick_max',max(abs(dv)),'force_rms',norm(Fstart)/sqrt(m), ...
    'force_max',max(abs(Fstart)),'force_projection',real(parent_direction'*Fstart)/m, ...
    'force_change_projection',real(parent_direction'*(Fstart-Fparent))/m);
assert(abs(initial.kick_RMS-options.kick)<1e-12, ...
    'LB:HexFlowKick','The start does not have the configured RMS kick.');
N=LB.N; L=LB.L; par=LB.par; cfg=options; x=phi; v=parent_direction; init=initial; %#ok<NASGU>
save(initial_file,'x','v','init','cfg','-v7.3');
write_json(fullfile(output_dir,[tag '_initial.json']),initial);
fid=fopen(history_file,'w'); assert(fid>0);
fprintf(fid,'updates,force_rms,force_max,energy_if_recorded,elapsed_seconds\n');
history=struct('updates',{},'updates_lower_bound',{},'updates_upper_bound',{}, ...
    'force_rms',{},'force_max',{},'energy',{},'energy_if_recorded',{},'elapsed_seconds',{});
run_clock=tic; status='budget_exhausted'; err=''; updates=0;
residual=NaN; force_max=NaN; energy=NaN; zero=zeros(m,1);
try
    for n=0:options.maxiter
        updates=n;
        f=lb_negative_gradient(phi); residual=norm(f)/sqrt(m); force_max=max(abs(f));
        if ~all(isfinite(phi)) || ~isfinite(residual)
            error('LB:Nonfinite','Nonfinite GF state');
        end
        if mod(n,options.energy_every)==0 || residual<options.tol || n==options.maxiter
            energy=lb_energy(phi); sample_energy=energy;
        else
            sample_energy=NaN;
        end
        sample_seconds=toc(run_clock);
        fprintf(fid,'%d,%.17g,%.17g,%.17g,%.9f\n',n,residual,force_max,sample_energy,sample_seconds);
        history(end+1)=struct('updates',n,'updates_lower_bound',n,'updates_upper_bound',n, ...
            'force_rms',residual,'force_max',force_max,'energy',sample_energy, ...
            'energy_if_recorded',sample_energy,'elapsed_seconds',sample_seconds); %#ok<AGROW>
        if mod(n,options.flush_every)==0
            fclose(fid); fid=fopen(history_file,'a'); assert(fid>0);
            fprintf('GF sign=%+d updates=%d force_RMS=%.6e force_max=%.6e energy=%.12g elapsed=%.2f\n', ...
                branch_sign,n,residual,force_max,energy,toc(run_clock));
        end
        if residual<options.tol, status='force_converged'; break; end
        if n==options.maxiter, break; end
        phi=position_update(phi,zero,options.dt);
        updates=n+1;
        if mod(n+1,options.checkpoint_every)==0, write_checkpoint(toc(run_clock)); end
    end
catch exception
    status='failed'; err=getReport(exception,'extended','hyperlinks','off');
end
fclose(fid); solver_seconds=toc(run_clock); audit_clock=tic;
try
    f=lb_negative_gradient(phi); residual=norm(f)/sqrt(m);
    force_max=max(abs(f)); energy=lb_energy(phi);
catch exception
    status='failed'; err=[err newline getReport(exception,'extended','hyperlinks','off')];
end
endpoint_audit_seconds=toc(audit_clock); elapsed_seconds=toc(run_clock);
write_checkpoint(elapsed_seconds);
flow=struct('status',status,'converged',strcmp(status,'force_converged'), ...
    'updates',updates,'updates_lower_bound',updates,'updates_upper_bound',updates, ...
    'update_count_exact',true,'force_rms',residual,'force_max',force_max,'energy',energy, ...
    'elapsed_seconds',elapsed_seconds,'solver_seconds',solver_seconds, ...
    'endpoint_audit_seconds',endpoint_audit_seconds,'history',history, ...
    'parent_phi_sha256',parent_hash,'direction_sha256',direction_hash, ...
    'initial_phi_sha256',initial_hash,'submitted_initial_phi_sha256',submitted_initial_hash, ...
    'branch_sign',branch_sign,'options',options,'initial',initial, ...
    'checkpoint_file',checkpoint_file,'history_file',history_file,'endpoint_file',endpoint_file, ...
    'resumed',false,'minimum_hessian_verified',false,'error',err);
x=phi; branch=flow; branch=rmfield(branch,'history'); %#ok<NASGU>
save(endpoint_file,'x','phi','v','cfg','branch','flow','-v7.3');
write_json(fullfile(output_dir,[tag '_endpoint.json']),branch);
fprintf('GF_END sign=%+d status=%s n=%d force=%.6e f=%.12g seconds=%.2f\n', ...
    branch_sign,status,updates,residual,energy,elapsed_seconds);

    function write_checkpoint(seconds)
        checkpoint=struct('phi',phi,'updates',updates,'lower',updates,'upper',updates, ...
            'history',history,'elapsed_seconds',seconds,'parent_hash',parent_hash, ...
            'direction_hash',direction_hash,'initial_phi_hash',initial_hash, ...
            'branch_sign',branch_sign,'gf_options',options,'N',N,'L',L,'par',par);
        tmp=fullfile(output_dir,[tag '_checkpoint.tmp.mat']);
        save(tmp,'-struct','checkpoint','-v7.3');
        movefile(tmp,checkpoint_file,'f');
    end
end

function phi=position_update(phi,reflected,b)
% Preserve the np_model('position') expression and its floating-point order.
global LB
phi=reshape(phi,LB.N,LB.N,LB.N); phi=phi-mean(phi(:));
reflected=reshape(reflected,size(phi)); p=LB.par;
low=p.tau*phi-(p.gamma/2)*phi.^2+phi.^3/6;
z=real(ifftn(fftn(phi+b*(-low-reflected))./(LB.ones_hat+b*LB.Lhigh_hat)));
z=z-mean(z(:)); phi=z(:);
end

function options=with_defaults(options,defaults)
names=fieldnames(defaults);
for j=1:numel(names)
    if ~isfield(options,names{j}) || isempty(options.(names{j}))
        options.(names{j})=defaults.(names{j});
    end
end
end

function write_json(file,value)
fid=fopen(file,'w'); assert(fid>0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
