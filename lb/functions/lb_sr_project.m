function Y=lb_sr_project(X,Q)
% Mean-zero plus low-rank RMS-orthogonal constraint, never a dense projector.
Y=X-mean(X,1);
if ~isempty(Q); Y=Y-Q*(Q'*Y/size(X,1)); end
end
