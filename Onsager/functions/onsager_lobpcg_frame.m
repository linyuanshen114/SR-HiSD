function [V,Vp] = onsager_lobpcg_frame(next,oldV,previous_V,nextU,p,target)
%ONSAGER_LOBPCG_FRAME Original three-block Ritz update of the ascent frame.
HV = onsager_hessian_action(next,oldV,nextU,p);
alpha = diag(oldV'*HV*p.w)';
residual = onsager_project((HV-oldV.*alpha).*exp(next),nextU,p.w);
trial = [oldV,residual,onsager_project(previous_V,nextU,p.w)];
trial = onsager_orthogonalize(trial,p.w);
assert(size(trial,2)>=target,'LOBPCG trial subspace lost rank.');
Htrial = onsager_hessian_action(next,trial,nextU,p);
small = trial'*Htrial*p.w;
[eta,d] = eig((small+small')/2,'vector');
[~,order] = sort(real(d));
V = onsager_orthogonalize(onsager_project(trial*real(eta(:,order(1:target))),nextU,p.w),p.w);
Vp = oldV;
end
