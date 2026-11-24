function certificate = morse_certificate(q,qmin,c)
%MORSE_CERTIFICATE Independently certify an endpoint in the rigid-body quotient.
% Index is counted on the SVD complement of two translations and one rotation.
certificate=struct('E',NaN,'gn',Inf,'rms',Inf,'index',-1, ...
    'components',0,'converged',false,'success',false,'mother',false,'mother_error',Inf);
if any(~isfinite(q)), return; end
[E,g,H]=morse_energy_gradient_hessian(q,c.a,c.r0);
N=numel(q)/2;
positions=reshape(q,2,N)'; centered=positions-mean(positions,1);
rotation=reshape([-centered(:,2),centered(:,1)]',[],1);
translations=repmat(eye(2),N,1);
Q=null([translations,rotation]');
Hphysical=Q'*H*Q;
certificate.E=E; certificate.gn=norm(g); certificate.rms=norm(g)/sqrt(numel(q));
certificate.index=nnz(eig((Hphysical+Hphysical')/2)<-c.eig_neg_tol);
% Connectivity is a finite graph test at the same distance cutoff as Table 1.
dx=positions(:,1)-positions(:,1)'; dy=positions(:,2)-positions(:,2)';
adjacency=dx.^2+dy.^2<=c.cluster_cutoff^2;
certificate.components=max(conncomp(graph(adjacency,'omitselfloops')));
certificate.converged=isfinite(E) && certificate.rms<=c.res_tol && certificate.gn<=c.grad_tol;
certificate.success=certificate.converged && certificate.components==1 && certificate.index==1;
% A mother return must align the same particle labels by a proper rotation.
if certificate.converged && certificate.index==0
    reference=reshape(qmin,2,N)'; reference=reference-mean(reference,1);
    [U,~,W]=svd(centered'*reference);
    R=U*diag([1,sign(det(U*W'))])*W';
    certificate.mother_error=norm(centered*R-reference,'fro');
    certificate.mother=certificate.components==1 && certificate.mother_error<=1e-6;
end
end
