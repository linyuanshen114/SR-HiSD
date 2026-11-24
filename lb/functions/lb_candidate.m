function c=lb_candidate(phi,solver,target,index_tolerance,force_tolerance,eigenpair_count)
%LB_CANDIDATE Independent coarse full/reduced Hessian diagnostics.
% These 100-step LOBPCG counts do not constitute strict index certificates.
% HEX direct SR uses eight pairs; native LBSD and disorder use ten.
if nargin<6, eigenpair_count=10; end
[V_reduced,lambda_reduced]=lb_reduced_eigenpairs(phi,eigenpair_count);
[V_full,lambda_full]=lb_full_eigenpairs(phi,eigenpair_count);
[lambda_full,order]=sort(real(lambda_full(:))); V_full=V_full(:,order);
force=lb_negative_gradient(phi); residual=sqrt((force'*force)/numel(phi));
index=sum(lambda_reduced < -index_tolerance); full_index=sum(lambda_full < -index_tolerance);
if isempty(target), target=index; end
eigen_residual_full=zeros(numel(lambda_full),1);
for j=1:numel(lambda_full)
    residual_vector=lb_hessian_action(phi,V_full(:,j))-lambda_full(j)*V_full(:,j);
    eigen_residual_full(j)=sqrt((residual_vector'*residual_vector)/numel(phi));
end
c=struct('phi',phi,'out',solver,'V_reduced',V_reduced,'lambda',lambda_reduced(:), ...
    'V_full',V_full,'lambda_full',lambda_full,'index',index,'full_index',full_index, ...
    'eigen_residual_full',eigen_residual_full, ...
    'full_zero_count',sum(abs(lambda_full)<=index_tolerance), ...
    'spectra_finite',all(isfinite(lambda_reduced)) && all(isfinite(lambda_full)), ...
    'residual',residual,'energy',lb_energy(phi),'target',target, ...
    'valid',all(isfinite(phi)) && residual<=1.05*force_tolerance && ...
        index==target,'strict_index_certified',false);
fprintf('Candidate: index=%d full=%d target=%d RMS=%.3e E=%.12g valid=%d\n', ...
    index,full_index,target,residual,c.energy,c.valid);
end
