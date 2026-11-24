function w = lb_precondition(w)
% lb_precondition(w)
% Simple spectral preconditioner approximating inverse of the stiff linear operator.

global LB
N = LB.N;
par = LB.par;

w3 = reshape(w, N,N,N);
w3 = w3 - mean(w3(:));

What = fftn(w3);

% Approx inverse for (Lhigh + shift)
shift = abs(par.tau) + 1.0; % mild positive shift
Pinv_hat = 1 ./ (LB.Lhigh_hat + shift);

w3p = real(ifftn(Pinv_hat .* What));
w3p = w3p - mean(w3p(:));
w = w3p(:);
end
