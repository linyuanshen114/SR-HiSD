function [x, output] = lb_lbsd(F, H, x, options)
% lb_lbsd(F,H,x,options)
% Semi-implicit HiSD / 1-LBSD style iteration for 3D LB model (CAM/FFT)
%
% F(x): negative gradient (vector)
% H(x,v): Hessian-vector product
% x: current state (vector, length N^3)
%
% options:
%   innpfunc, precd, dt, maxiter, epsf, betat, betau
%   k, V (initial directions)
%   outputp, outputd
%   drawfun (optional) : @(x,it,res,alph) ...

global LB
N = LB.N;
solver_tic = tic;

% inner product
if isfield(options,'innpfunc') && ~isempty(options.innpfunc)
    innp = @(X,Y) options.innpfunc(X,Y);
    nor  = @(v) sqrt(innp(v,v));
else
    innp = @(X,Y) (X'*Y);
    nor  = @(v) norm(v);
end

% preconditioner
if isfield(options,'precd') && ~isempty(options.precd)
    P = options.precd;
else
    P = @(w) w;
end

dt = options.dt;
epsf = options.epsf;
maxiter = options.maxiter;
k = options.k;
V = options.V;
% Shared with SR-HiSD: discard only numerically degenerate trial blocks.
block_tol = 1e-8;

% Orthonormalize V
for i=1:k
    V(:,i) = V(:,i) - V(:,1:i-1)*innp(V(:,1:i-1), V(:,i));
    V(:,i) = V(:,i) / max(nor(V(:,i)), 1e-14);
end

% BB helpers
bta = dt;
xp = x;
gp = zeros(size(x));

alph = zeros(k,1);
U = zeros(numel(x),k);
W = zeros(numel(x),k);
Y = zeros(numel(x),k);
Vp = zeros(numel(x),k);
Up = zeros(numel(x),k);

for it = 1:maxiter
    % Force and reflected force
    f = F(x);
    tmpp = 2 * V * innp(V, f);     % 2*Proj_V(f)
    g = f - tmpp;                  % HiSD reflected force (descend stable / ascend along V)

    if it > 1
        Dx = x - xp;
        Dg = g - gp;
        denom = innp(Dg,Dg);
        if abs(denom) > 1e-20
            bta = abs(innp(Dx,Dg) / denom);
        else
            bta = dt;
        end
        % clamp
        bta = min(bta, options.betat*dt);
        bta = max(bta, options.betau*dt);
    end
    xp = x;
    gp = g;

    % -------- semi-implicit update on the stiff linear part (Lhigh) --------
    phi3 = reshape(x, N,N,N);
    phi3 = phi3 - mean(phi3(:));

    par = LB.par;
    grad_low = par.tau*phi3 - (par.gamma/2)*phi3.^2 + (1/6)*phi3.^3;

    tmpp3 = reshape(tmpp, N,N,N);

    rhs = phi3 + bta * ( -grad_low - tmpp3 ); % explicit low/nonlinear + reflected part
    rhs_hat = fftn(rhs);

    denom_hat = LB.ones_hat + bta * LB.Lhigh_hat;
    phi_new = real(ifftn(rhs_hat ./ denom_hat));
    phi_new = phi_new - mean(phi_new(:));
    x = phi_new(:);
    % -------- update unstable directions V via LOBPCG-like step --------
    emp = true(3*k,1);

    for i=1:k
        U(:,i) = H(x, V(:,i));
        alph(i) = innp(V(:,i), U(:,i));

        W(:,i) = U(:,i) - alph(i) * V(:,i);
        W(:,i) = P(W(:,i));
        W(:,i) = W(:,i) - V * innp(V, W(:,i));
        W(:,i) = W(:,i) - W(:,1:i-1) * innp(W(:,1:i-1), W(:,i));

        nrmW = nor(W(:,i));
        if nrmW > block_tol
            W(:,i) = W(:,i)/nrmW;
            Y(:,i) = H(x, W(:,i));
        else
            W(:,i) = 0; Y(:,i) = 0;
            emp(k+i) = false;
        end
    end

    for i=1:k
        zeta = innp([V W Vp(:,1:i-1)], Vp(:,i));
        Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * zeta;
        nrmVpi = nor(Vp(:,i));
        if nrmVpi > block_tol
            Vp(:,i) = Vp(:,i)/nrmVpi;
            Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * innp([V W Vp(:,1:i-1)], Vp(:,i));
            Vp(:,i) = Vp(:,i) / max(nor(Vp(:,i)), 1e-14);
            Up(:,i) = H(x, Vp(:,i));
        else
            Vp(:,i) = 0; Up(:,i) = 0;
            emp(2*k+i) = false;
        end
    end

    UU = [V, W, Vp]; UU = UU(:, emp);
    YY = [U, Y, Up]; YY = YY(:, emp);

    Pn = innp(UU, YY); Pn = (Pn + Pn')/2;
    [eta, eigs_vec] = eig(Pn, 'vector');
    [eigs_vec, ord] = sort(real(eigs_vec), 'ascend');
    eta = eta(:, ord);

    Vp = V;
    V = UU * real(eta(:,1:k));
    alph = real(eigs_vec(1:k));

    % normalize V again
    for i=1:k
        V(:,i) = V(:,i) - V(:,1:i-1)*innp(V(:,1:i-1), V(:,i));
        V(:,i) = V(:,i) / max(nor(V(:,i)), 1e-14);
    end

    % residual
    fnew = F(x);
    res = sqrt(innp(fnew,fnew));

    if mod(it, options.outputp) == 0
        alphstr = sprintf('%.6f ', alph);
        fprintf('LBSD it=%d  res=%.3e  alph=%s  bta=%.3g\n', it, res, alphstr, bta);
    end

    if isfield(options,'drawfun') && ~isempty(options.drawfun) && mod(it, options.outputd) == 0
        options.drawfun(x, it, res, alph);
    end

    if res < epsf
        break;
    end
end

output = struct('x', x, 'V', V, 'it', it, 'alph', alph, 'res', res, ...
    'lobpcg_block_tol', block_tol, 'solve_seconds', toc(solver_tic));
end
