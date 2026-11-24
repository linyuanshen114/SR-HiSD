function U = onsager_constraint_basis(logrho,p)
c = p.L'\(p.L\(p.Y'*(logrho*p.w)));
rho = exp(logrho);
R = zeros(p.N,3);
for j = 1:3
    R(:,j) = rho.*((p.Ym{j}*c-p.Yp{j}*c)/(2*p.epsrot));
    R(:,j) = R(:,j)-mean(R(:,j));
end
A = [ones(p.N,1),R]; U = zeros(p.N,0);
for j = 1:size(A,2)
    v = onsager_project(A(:,j),U,p.w);
    nv = onsager_wnorm(v,p.w);
    if nv > p.sr_tol*max(1,onsager_wnorm(A(:,j),p.w)), U(:,end+1) = v/nv; end
end
end

