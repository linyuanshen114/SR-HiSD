function x = onsager_wnorm(v,w)
x = sqrt(max(0,v(:)'*v(:)*w));
end

