function [route,endpoints]=lb_reference_hex_lbsd(data_dir)
%LB_REFERENCE_HEX_LBSD Load the historical 1490068 sequential route and GF.
% Four accepted original-space stages are paired only with GF job 1490193.
base=fullfile(data_dir,'reference','bcc2hex_lbsd');
search=fullfile(base,'search'); gf=fullfile(base,'gf');
status=jsondecode(fileread(fullfile(search,'status.json')));
assert(strcmp(status.status,'original_space_index1_numerically_resolved'));
labels={'up_4','down_4to3_kick_1000_sign_+1', ...
    'down_3to2_kick_1000_sign_+1','down_2to1_kick_1000_sign_-1'};
indices=[4 3 2 1]; stages=cell(1,4);
for j=1:4
    file=fullfile(search,[labels{j} '_candidate.mat']);
    phi=double(h5read(file,'/c/phi')); phi=phi(:);
    index=double(h5read(file,'/c/actual_index'));
    accepted=logical(h5read(file,'/c/accepted'));
    residual=double(h5read(file,'/c/force_rms'));
    G=load(fullfile(search,[labels{j} '_original_gate.mat']),'G'); G=G.G;
    assert(accepted && index==indices(j) && residual<=1e-7 && ...
        all(isfinite(phi)) && G.verified && G.index==indices(j) && G.zero_tol==1e-7, ...
        'LBReview:SequentialGate','A saved sequential stage failed its original-space numerical index gate.');
    assert(numel(G.starts)==2);
    for k=1:2
        evidence=G.starts{k};
        assert(evidence.verified && evidence.index==indices(j) && evidence.zero_tol==1e-7 && ...
            nnz(evidence.signs==-1)==indices(j), ...
            'LBReview:SequentialGate','Warm/cold original-space numerical indices disagree.');
    end
    stages{j}=phi;
end
assert(isequal(indices,[4 3 2 1]));
selected_phi=stages{end};
selected_file=fullfile(search,[labels{end} '_candidate.mat']);
spectrum_file=fullfile(search,[labels{end} '_original_start1.mat']);
S=load(fullfile(gf,'source.mat'),'source','v','cfg'); source=S.source; v=S.v(:); cfg=S.cfg;
assert(strcmp(source.candidate_sha256,file_sha(selected_file)) && ...
    strcmp(source.spectrum_sha256,file_sha(spectrum_file)), ...
    'LBReview:GFSource','Saved GF has different candidate/spectrum file bytes.');
parent_hash=lb_numeric_sha(selected_phi); direction_hash=lb_numeric_sha(v);
assert(strcmp(source.phi_sha256,parent_hash) && strcmp(source.direction_sha256,direction_hash));
assert(source.index==1 && source.lambda<-1e-7 && source.eigen_residual<=1e-8 && source.force_rms<=1e-7);
Q=load(fullfile(data_dir,'input','bcc2hex_lbsd_reference.mat'),'input'); input=Q.input;
assert(input.N==128 && abs(input.L-24*sqrt(2)*pi)<1e-12 && ...
    input.par.tau==-.3 && input.par.gamma==.8 && input.par.xi==1 && input.par.q0==1);
assert(strcmp(source.parent_sha256,lb_numeric_sha(input.phi_parent)) && ...
    strcmp(source.hex_reference_sha256,lb_numeric_sha(input.hex_reference)));
assert(cfg.kick==1e-5 && cfg.dt==1 && cfg.tol==1e-9 && ~cfg.nesterov);
global LB; prior_model=LB;
model_cleanup=onCleanup(@() restore_model(prior_model));
lb_initialize(input.N,input.L,input.par);
force=lb_negative_gradient(selected_phi); Hv=lb_hessian_action(selected_phi,v);
lambda=real(v'*Hv)/numel(v); raw_residual=norm(Hv-lambda*v)/sqrt(numel(v));
assert(norm(force)/sqrt(numel(v))<=1e-7 && abs(lambda-source.lambda)<1e-8 && raw_residual<=1e-8, ...
    'LBReview:GFDirection','Saved negative mode failed a fresh force/raw-HVP audit.');
endpoints=struct([]); signs=[1 -1];
for j=1:2
    branch_dir=fullfile(gf,sprintf('sign_%+d',signs(j)));
    D=load(fullfile(branch_dir,'endpoint.mat'),'x','v','branch','cfg');
    I=load(fullfile(branch_dir,'initial.mat'),'x','v','init','cfg');
    assert(isequaln(D.v(:),v) && isequaln(I.v(:),v) && isequaln(D.cfg,cfg) && isequaln(I.cfg,cfg));
    branch=D.branch; phase=branch.phase_evidence;
    assert(branch.sign==signs(j) && I.init.sign==signs(j) && ...
        strcmp(branch.status,'force_converged') && branch.force_rms<1e-9);
    expected=selected_phi+signs(j)*cfg.kick*v; expected=expected-mean(expected);
    assert(norm(I.x(:)-expected)/sqrt(numel(v))<1e-13 && max(abs(I.x(:)-expected))<1e-12, ...
        'LBReview:GFInitialState','Saved GF initial field is not its stated parent/sign/direction.');
    if strcmp(phase.phase,'HEX')
        assert(phase.hex_reference_match.verified && phase.primary_rank==2 && ...
            phase.six_peak_power>.99 && phase.six_peak_relative_spread<1e-5 && ...
            phase.offplane_power<1e-10 && abs(phase.triad_phase_closure)<1e-6);
    elseif strcmp(phase.phase,'BCC')
        assert(phase.parent_reference_match.verified);
    else
        error('LBReview:EndpointPhase','The sequential GF phase is unresolved.');
    end
    hex_reference=struct('phi',input.hex_reference);
    endpoint=lb_branch_info(D.x(:),input.phi_parent(:),selected_phi,signs(j),hex_reference);
    assert(endpoint.converged && strcmp(endpoint.class,phase.phase), ...
        'LBReview:EndpointPhase','A fresh endpoint phase audit disagrees with the saved phase evidence.');
    endpoint.source_phase_evidence=phase;
    endpoint.flow=struct('branch_sign',signs(j),'parent_phi_sha256',parent_hash, ...
        'direction_sha256',direction_hash,'initial_phi_sha256',lb_numeric_sha(I.x), ...
        'converged',true,'iterations',branch.updates,'iterations_exact',true, ...
        'elapsed_seconds',branch.elapsed_seconds,'solver_seconds',branch.solver_seconds,'options',cfg);
    if isempty(endpoints), endpoints=endpoint; else, endpoints(j)=endpoint; end
end
seed=struct('N',input.N,'L',input.L,'par',input.par,'kick',.01, ...
    'phi_bcc',input.phi_parent(:),'bb',struct('epsf',1e-7));
selected=struct('phi',selected_phi,'residual',source.force_rms,'index',1,'valid',true, ...
    'negative_mode_verified',true,'independent_index1_verified',true,'strict_index_certified',false);
route=lb_make_route('bcc2hex','lbsd4_down',seed,selected,stages,endpoints,false,indices);
assert(route.complete && route.plot_ready && nnz(strcmp({endpoints.class},'HEX'))==1 && ...
    nnz(strcmp({endpoints.class},'BCC'))==1);
route.fresh=false;
clear model_cleanup
end

function value=file_sha(file)
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(file,'rb'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
while true
    bytes=fread(fid,4*1024*1024,'*uint8');
    if isempty(bytes), break; end
    md.update(bytes);
end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
clear cleanup
end

function restore_model(model)
global LB; LB=model;
end
