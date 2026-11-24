function V = onsager_canonical_signs(V)
for j = 1:size(V,2)
    [~,i] = max(abs(V(:,j)));
    if V(i,j)<0, V(:,j)=-V(:,j); end
end
end

