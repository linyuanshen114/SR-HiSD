function routes=lb_reference_routes(data_dir,case_id,kicks)
%LB_REFERENCE_ROUTES Load the four selected saved fields and actual endpoints.
% Cached SR fields use only independently completed GF from the exact parent.
if nargin<3, kicks=.02; end
assert(all(abs(kicks-.02)<1e-12),'LBReview:Kick','The selected disorder case uses kick 0.02.');
switch case_id
    case 'bcc2hex'
        [lb,endpoints]=lb_reference_hex_lbsd(data_dir);
    case 'bcc2disorder'
        file='archived/bcc2disorder/kick_0020/lbsd4_to1/route_final.mat';
        S=lb_read_view(data_dir,file,{'N','L','par','bb','kick','selected','phi_bcc', ...
            'phi_plus','phi_minus','endpoint_plus','endpoint_minus'});
        E={S.endpoint_plus,S.endpoint_minus}; P={S.phi_plus,S.phi_minus}; signs=[1 -1]; endpoints=struct([]);
        for j=1:2
            endpoint=struct('phi',P{j},'sign',signs(j),'class',E{j}.class, ...
                'residual',E{j}.residual,'converged',E{j}.residual<1e-9, ...
                'parent_phi_sha256',lb_numeric_sha(S.selected.phi),'phi_sha256',lb_numeric_sha(P{j}));
            if isempty(endpoints), endpoints=endpoint; else, endpoints(j)=endpoint; end
        end
        seed=struct('N',S.N,'L',S.L,'par',S.par,'kick',S.kick,'phi_bcc',S.phi_bcc,'bb',S.bb);
        selected=selection(S.selected.phi,S.selected.residual,S.selected.morse);
        lb=lb_make_route(case_id,'lbsd4_to1',seed,selected,{S.selected.phi},endpoints,false,S.selected.morse);
    otherwise
        error('LBReview:Case','Unknown case.');
end
assert(lb.complete && all(strcmp({endpoints.parent_phi_sha256},lb.parent_phi_sha256)), ...
    'LBReview:SavedEndpoints','Saved LBSD endpoints must be complete and share the exact parent.');
lb.fresh=false;

%% Each cached SR field is loaded from its own exact input and final state.
name=[case_id '_srhisd'];
I=load(fullfile(data_dir,'input',[name '.mat']),'input'); input=I.input;
F=load(fullfile(data_dir,'reference',name,'final.mat'),'state');
meta=jsondecode(fileread(fullfile(data_dir,'reference',name,'metadata.json')));
summary=jsondecode(fileread(fullfile(data_dir,'reference',name,'summary.json')));
assert(strcmp(lb_numeric_sha(input.phi_init),meta.source_hashes.phi_init));
assert(strcmp(lb_numeric_sha(input.V_seed),meta.source_hashes.V_seed));
seed=struct('N',input.N,'L',input.L,'par',input.par,'kick',meta.configuration.kick, ...
    'phi_bcc',input.phi_parent,'bb',struct('epsf',meta.configuration.tol));
selected=selection(F.state.phi,summary.force_rms,NaN);
if strcmp(case_id,'bcc2hex')
    proof=jsondecode(fileread(fullfile(data_dir,'reference',name,'certification','verdict.json')));
    V=load(fullfile(data_dir,'reference',name,'verified_negative.mat'),'source_phi_sha256');
    assert(strcmp(char(V.source_phi_sha256),lb_numeric_sha(F.state.phi)));
    selected.independent_index1_verified=proof.verified_index1;
    selected.negative_mode_verified=true;
end
%% Load a completed same-parent GF pair when it has been supplied.
sr_endpoints=struct([]);
gf_dir=fullfile(data_dir,'reference',name,'gf');
completion_file=fullfile(gf_dir,'completion.json');
if isfile(completion_file)
    completion=jsondecode(fileread(completion_file));
    parent_hash=lb_numeric_sha(F.state.phi);
    assert(strcmp(completion.status,'PASS') && strcmp(completion.parent_phi_sha256,parent_hash), ...
        'LBReview:GFParent','Saved GF belongs to another SR parent.');
    assert(completion.negative_mode_verified && completion.saddle_search_repeated==false);
    assert(completion.N==seed.N && abs(completion.L-seed.L)<1e-12);
    for parameter={'tau','gamma','xi','q0'}
        assert(completion.par.(parameter{1})==seed.par.(parameter{1}), ...
            'LBReview:GFModel','Saved GF uses different model parameters.');
    end
    selected.negative_mode_verified=true;
    names={'plus','minus'}; signs=[1 -1];
    for j=1:2
        D=load(fullfile(gf_dir,names{j},'endpoint.mat'),'endpoint');
        endpoint=D.endpoint;
        assert(endpoint.sign==signs(j) && endpoint.converged && endpoint.residual<1e-9);
        assert(strcmp(endpoint.parent_phi_sha256,parent_hash) && ...
            strcmp(endpoint.flow.parent_phi_sha256,parent_hash), ...
            'LBReview:GFParent','The GF endpoint has a different exact parent.');
        assert(strcmp(endpoint.phi_sha256,lb_numeric_sha(endpoint.phi)) && ...
            strcmp(endpoint.flow.direction_sha256,completion.direction_sha256));
        assert(endpoint.flow.branch_sign==signs(j) && endpoint.flow.converged);
        assert(endpoint.flow.options.dt==1 && endpoint.flow.options.tol==1e-9);
        if strcmp(endpoint.class,'HEX')
            assert(isfield(endpoint,'phase_evidence') && endpoint.phase_evidence.hex_phase_certified, ...
                'LBReview:HEXPhase','HEX requires independent phase evidence.');
        end
        if isempty(sr_endpoints), sr_endpoints=endpoint; else, sr_endpoints(j)=endpoint; end
    end
end
sr=lb_make_route(case_id,'srhisd1',seed,selected,{F.state.phi},sr_endpoints,false,1);
sr.fresh=false;
assert(sr.plot_ready,'LBReview:SavedSR','The saved SR field failed the force tolerance.');
if isfile(completion_file)
    assert(sr.complete,'LBReview:SavedConnection','The supplied SR GF pair does not establish the expected phase connection.');
end
routes=[lb sr];
end

function s=selection(phi,residual,index)
s=struct('phi',phi,'residual',residual,'valid',isfinite(residual) && residual<=1e-7, ...
    'index',index,'strict_index_certified',false,'negative_mode_verified',false, ...
    'independent_index1_verified',false);
end
