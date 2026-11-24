function Q=lb_sr_orth(X,project,tol)
if nargin<3; tol=1e-12; end
Q=zeros(size(X,1),0); X=project(X); m=size(X,1);
for j=1:size(X,2)
    z=X(:,j);
    for pass=1:2; z=z-Q*(Q'*z/m); end
    z=project(z); n=norm(z)/sqrt(m);
    if n>tol; Q(:,end+1)=z/n; end %#ok<AGROW>
end
end
