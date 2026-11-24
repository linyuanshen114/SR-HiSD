function [z,error] = onsager_entropy_solve(beta,rhs,z)
% Monotone scalar equations exp(z)+beta*z=rhs, one equation per node.
threshold = 2e-14*max(1,max(abs(rhs)));
for j = 1:60
    rho = exp(z);
    r = rho+beta*z-rhs;
    if max(abs(r)) <= threshold, break; end
    dz = r./(rho+beta);
    z = z-max(-4,min(4,dz));
end
error = max(abs(exp(z)+beta*z-rhs));
assert(error<=5*threshold,'Scalar entropy solve failed.');
end

