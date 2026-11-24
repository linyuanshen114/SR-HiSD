function C=lb_classify_zero_band(lambda,residual,ztol)
%LB_CLASSIFY_ZERO_BAND Classify raw-Hessian Ritz residual intervals.
% Signs are -1 negative, 0 in the zero band, +1 positive, 9 unresolved.
lambda=lambda(:); residual=residual(:); assert(numel(lambda)==numel(residual));
assert(all(isfinite(lambda)) && all(isfinite(residual)) && all(residual>=0));
lo=lambda-residual; hi=lambda+residual; s=9*ones(size(lambda));
s(hi < -ztol)=-1; s(lo > ztol)=1; s(lo>=-ztol & hi<=ztol)=0;
pos=find(s==1); C=struct('verified',false,'index',NaN,'prefix',0,'signs',s,'lo',lo,'hi',hi,'zero_tol',ztol);
if numel(pos)>=2
    p=pos(2); C.prefix=p;
    if all(s(1:p)~=9) && all(residual(1:p)<=1e-8)
        C.verified=true; C.index=nnz(s(1:p)==-1);
    end
end
end
