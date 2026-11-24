function E = lb_energy(phi, phihat)
% lb_energy(phi)
% Mean LB free-energy density on the periodic cube [0,L]^3.

global LB
N = LB.N; L = LB.L;
par = LB.par;

phi3 = reshape(phi, N,N,N);
phi3 = phi3 - mean(phi3(:));

if nargin < 2 || isempty(phihat); phihat = fftn(phi3); end

% (Δ+q0^2)phi in Fourier: (q0^2 - |q|^2) * phihat
lap1 = real(ifftn((par.q0^2 - LB.K2) .* phihat));
term_grad = 0.5 * par.xi^2 * mean(lap1(:).^2);

term_quad = 0.5 * par.tau * mean(phi3(:).^2);
term_cub  = -(par.gamma/6) * mean(phi3(:).^3);
term_qua  = (1/24) * mean(phi3(:).^4);

% Return the spatial mean, matching the staged solver's energy convention.
E = (term_grad + term_quad + term_cub + term_qua);

end
