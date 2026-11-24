function V = onsager_orthogonalize(A,w)
V = zeros(size(A,1),0);
for j = 1:size(A,2)
    v = A(:,j);
    for pass = 1:2, v = onsager_project(v,V,w); end
    nv = onsager_wnorm(v,w);
    if nv>1e-10, V(:,end+1) = v/nv; end
end
end

