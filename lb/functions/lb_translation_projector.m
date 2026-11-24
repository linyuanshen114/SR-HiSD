function proj = lb_translation_projector(phi, cfg, phihat)
% lb_translation_projector(phi, cfg)
% Build operator-form orthogonal projectors onto:
%   - tangent space of translation orbit: span{dphi/dx, dphi/dy, dphi/dz}
%     or a user-specified subset/linear combination via cfg.directions.
%   - normal space (quotient direction): I - P_tan
%
% proj.PT(v): project v to tangent space
% proj.PN(v): project v to normal space
% proj.Q    : L2-orthonormal tangent basis (d x r, r<=3)
% proj.rank : r
% proj.G    : Gram matrix (3x3)

if nargin < 2 || isempty(cfg), cfg = struct(); end
if ~isfield(cfg,'tol') || isempty(cfg.tol), cfg.tol = 1e-6; end

global LB
N = LB.N;
d = numel(phi);

phi3 = reshape(phi, N,N,N);
if nargin < 3 || isempty(phihat), phihat = fftn(phi3); end

% spectral derivatives: dphi/dx, dphi/dy, dphi/dz  (Fourier: multiply by i*k)
tx = real(ifftn(1i * LB.KX .* phihat));
ty = real(ifftn(1i * LB.KY .* phihat));
tz = real(ifftn(1i * LB.KZ .* phihat));

T = [tx(:), ty(:), tz(:)];
if isfield(cfg, 'directions') && ~isempty(cfg.directions)
    dirs = cfg.directions;
    if size(dirs, 1) ~= 3 && size(dirs, 2) == 3
        dirs = dirs.';
    end
    if size(dirs, 1) ~= 3
        error('cfg.directions must have size 3-by-r or r-by-3.');
    end
    T = T * dirs;
end

% Gram matrix under <u,v> = (u'*v)/d
G = (T.' * T) / d;
G = 0.5*(G + G.');

% rank-adaptive orthonormalization via eig(G) (3x3 only)
[V,Dm] = eig(G);
lam = real(diag(Dm));
[lam, idx] = sort(lam, 'descend');
V = V(:, idx);

lam0 = max(lam(1), 0);
keep = lam > cfg.tol * max(lam0, 1);
r = sum(keep);

if r == 0
    Q = zeros(d,0);
else
    Vk = V(:, keep);
    lk = lam(keep);
    % Q = T*Vk*diag(1/sqrt(lk))  => Q'*Q/d = I
    Q = T * Vk * diag(1 ./ sqrt(lk));
end

proj = struct();
proj.Q = Q;
proj.rank = r;
proj.G = G;
if isfield(cfg, 'directions')
    proj.directions = cfg.directions;
end

proj.PT = @(v) apply_PT(v, Q, d);
proj.PN = @(v) (v - apply_PT(v, Q, d));

end

function y = apply_PT(v, Q, d)
% Apply tangent projection using L2-orthonormal Q: PT(v) = Q*(Q'*v/d)
if isempty(Q)
    y = zeros(size(v));
    return;
end

V = v;
isvec = isvector(V);
if isvec, V = V(:); end

y = Q * ((Q.' * V) / d);

if isvec, y = y(:); end
end
