function Y = onsager_spherical_matrix(xyz,lmax)
xyz = xyz./vecnorm(xyz,2,2);
ct = cos(acos(max(-1,min(1,xyz(:,3)))));
phi = atan2(xyz(:,2),xyz(:,1));
Y = zeros(size(xyz,1),(lmax+1)^2); col = 0;
for l = 0:lmax
    P = legendre(l,ct');
    for m = -l:l
        col = col+1; mm = abs(m);
        c = sqrt((2*l+1)/(4*pi)*exp(gammaln(l-mm+1)-gammaln(l+mm+1)));
        y = c*P(mm+1,:)'.*exp(1i*mm*phi);
        if m<0
            Y(:,col)=sqrt(2)*imag(y);
        elseif m==0
            Y(:,col)=real(y);
        else
            Y(:,col)=sqrt(2)*real(y);
        end
    end
end
end

