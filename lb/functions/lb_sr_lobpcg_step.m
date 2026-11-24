function [Vnew,Vprevious,info]=lb_sr_lobpcg_step(Hv,V,Vprevious,project,precondition,opts)
%LB_SR_LOBPCG_STEP Complete-projection generalized-Gram direction kernel.
% Both methods retain the SAME history and assembled B/HB projections.
% Local LB adaptation: RMS inner products, generalized-Gram whitening,
% two-pass GS, allowed-space trial vectors, and cached single-column HVPs.
% project must be a fixed linear orthogonal projector throughout this step.
% k=1 with three nonempty blocks: 21 callbacks / 25 scalar projected columns.
% project_seconds measures callback execution only, before finite checks.
% Rank-deficient historical columns are dropped, never replaced by QR filler.
if nargin<3 || isempty(Vprevious); Vprevious=[]; end
if nargin<4 || isempty(project); project=@(X) X; end
if nargin<5; precondition=[]; end
if nargin<6 || isempty(opts); opts=struct(); end
if ~isfield(opts,'block_tol'); opts.block_tol=1e-8; end
if ~isfield(opts,'gram_tol'); opts.gram_tol=max(100*eps,opts.block_tol^2); end
assert(isa(Hv,'function_handle') && isa(project,'function_handle'),'lb_sr_lobpcg_step:Callbacks','Invalid callbacks.');
assert(isempty(precondition) || isa(precondition,'function_handle'),'lb_sr_lobpcg_step:Callbacks','Invalid preconditioner.');
assert(isnumeric(V) && ismatrix(V) && ~isempty(V) && all(isfinite(V(:))),'lb_sr_lobpcg_step:InitialBlock','Invalid V.');
assert(isscalar(opts.block_tol) && isfinite(opts.block_tol) && opts.block_tol>=0 && ...
    isscalar(opts.gram_tol) && isfinite(opts.gram_tol) && opts.gram_tol>0,'lb_sr_lobpcg_step:Options','Invalid tolerances.');
[m,k]=size(V);
if isempty(Vprevious); Vprevious=zeros(m,0,'like',V); end
assert(isnumeric(Vprevious) && size(Vprevious,1)==m && all(isfinite(Vprevious(:))),'lb_sr_lobpcg_step:HistoryBlock','Invalid history.');
hvp_count=0; precondition_count=0; project_count=0; project_scalar_columns=0; project_seconds=0;
[V,kept]=orthogonalize(V,zeros(m,0,'like',V));
assert(numel(kept)==k,'lb_sr_lobpcg_step:InitialRank','Projected initial V lost rank.');
% Preserve the initial whole-history projection, in addition to the three
% per-column projections in orthogonalize_history below.
Vprevious=apply_project(Vprevious);
HV=hessian_action(V);
[rotation,old_values]=eig(hermitian(V'*HV/m),'vector');
[old_values,order]=sort(real(old_values),'ascend'); rotation=rotation(:,order);
V=V*rotation; HV=HV*rotation; Vold=V;
R=apply_project(HV-V*diag(old_values));
old_residual_rms=sqrt(sum(abs(R).^2,1)/m);
Wraw=zeros(size(R),'like',R);
for j=1:k
    if old_residual_rms(j)<=opts.block_tol; continue; end
    if isempty(precondition); Wraw(:,j)=R(:,j);
    else
        z=precondition(R(:,j)); precondition_count=precondition_count+1;
        assert(isequal(size(z),[m,1]) && all(isfinite(z(:))),'lb_sr_lobpcg_step:PreconditionerOutput','Invalid preconditioner output.');
        Wraw(:,j)=z;
    end
end
Wraw=apply_project(Wraw);
[W,~]=orthogonalize(Wraw,V);
[Pold,~]=orthogonalize_history(Vprevious,[V,W]);
% Deliberately KEEP these complete trial-block projections for both methods.
B=apply_project([V,W,Pold]);
HB=apply_project([HV,hessian_action(W),hessian_action(Pold)]);
G=hermitian(B'*B/m); T=hermitian(B'*HB/m);
[whitener,gram_rank]=gram_whitener(G);
assert(gram_rank>=k,'lb_sr_lobpcg_step:TrialRank','Trial space has fewer than k directions.');
[small_vectors,values]=eig(hermitian(whitener'*T*whitener),'vector');
[~,order]=sort(real(values),'ascend'); coefficients=whitener*small_vectors(:,order(1:k));
Vnew=apply_project(B*coefficients); HVnew=apply_project(HB*coefficients);
[normalizer,selected_rank]=gram_whitener(hermitian(Vnew'*Vnew/m));
assert(selected_rank==k,'lb_sr_lobpcg_step:SelectedRank','Selected block lost rank.');
Vnew=Vnew*normalizer; HVnew=HVnew*normalizer;
[rotation,values]=eig(hermitian(Vnew'*HVnew/m),'vector');
[values,order]=sort(real(values),'ascend'); rotation=rotation(:,order);
Vnew=Vnew*rotation; HVnew=HVnew*rotation; Vprevious=Vold;
residual=HVnew-Vnew*diag(values); leak=apply_project(Vnew)-Vnew;
info=struct('ritz_values',values(:),'residual_rms',sqrt(sum(abs(residual).^2,1)/m).', ...
    'orthogonality_error',norm(Vnew'*Vnew/m-eye(k),'fro'), ...
    'constraint_leak_rms',sqrt(sum(abs(leak).^2,1)/m).', ...
    'old_ritz_values',old_values(:),'old_residual_rms',old_residual_rms(:), ...
    'block_sizes',[k,size(W,2),size(Pold,2)],'trial_dimension',size(B,2),'gram_rank',gram_rank, ...
    'hvp_count',hvp_count,'precondition_count',precondition_count, ...
    'project_count',project_count,'project_scalar_columns',project_scalar_columns, ...
    'project_seconds',project_seconds,'status','ok');

    function Y=apply_project(X)
        if isempty(X); Y=X; return; end
        callback_clock=tic; Y=project(X); callback_seconds=toc(callback_clock);
        project_seconds=project_seconds+callback_seconds;
        project_count=project_count+1; project_scalar_columns=project_scalar_columns+size(X,2);
        assert(isequal(size(Y),size(X)) && all(isfinite(Y(:))),'lb_sr_lobpcg_step:ProjectorOutput','Invalid projector output.');
    end

    function Y=hessian_action(X)
        X=apply_project(X); Y=zeros(size(X),'like',X);
        for jj=1:size(X,2)
            z=Hv(X(:,jj)); hvp_count=hvp_count+1;
            assert(isequal(size(z),[m,1]) && all(isfinite(z(:))),'lb_sr_lobpcg_step:HessianOutput','Invalid HVP output.');
            Y(:,jj)=z;
        end
        Y=apply_project(Y);
    end

    function [Q,indices]=orthogonalize(C,against)
        C=apply_project(C); Q=zeros(m,0,'like',V); indices=zeros(1,0);
        for jj=1:size(C,2)
            w=C(:,jj);
            for pass=1:2
                for ii=1:size(against,2); w=w-against(:,ii)*(against(:,ii)'*w/m); end
                for ii=1:size(Q,2); w=w-Q(:,ii)*(Q(:,ii)'*w/m); end
            end
            w=apply_project(w); scale=norm(w)/sqrt(m);
            if scale>opts.block_tol; Q(:,end+1)=w/scale; indices(end+1)=jj; end %#ok<AGROW>
        end
    end

    function [Q,indices]=orthogonalize_history(C,against)
        % Initial whole-block P was already applied above. Each retained
        % column gets P, first GS/P/normalize, second GS/P/normalize.
        Q=zeros(m,0,'like',V); indices=zeros(1,0);
        for jj=1:size(C,2)
            w=apply_project(C(:,jj)); history_against=[against,Q];
            w=w-history_against*(history_against'*w/m);
            w=apply_project(w); scale=norm(w)/sqrt(m);
            if scale<=opts.block_tol; continue; end
            w=w/scale;
            w=w-history_against*(history_against'*w/m);
            w=apply_project(w); scale=norm(w)/sqrt(m);
            if scale>opts.block_tol; Q(:,end+1)=w/scale; indices(end+1)=jj; end %#ok<AGROW>
        end
    end

    function [Z,numerical_rank]=gram_whitener(small_G)
        [U,d]=eig(hermitian(small_G),'vector'); d=real(d);
        retain=d>opts.gram_tol*max(1,max(abs(d))); numerical_rank=nnz(retain);
        Z=U(:,retain)./sqrt(d(retain)).';
    end
end

function A=hermitian(A)
A=(A+A')/2;
end
