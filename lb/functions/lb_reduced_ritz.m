function [V, alpha, iter] = lb_reduced_ritz(Hvv,x,options)
%% reading parameters
n=length(x);
% innp can compute multi - one innp as a column vector
if isfield(options,'innpfunc')
    innp = @(x,y)options.innpfunc(x,y);
    nor = @(x)sqrt(innp(x,x));
else
    disp('HiOSD: standard inner product.');
    innp = @(x,y)(x'*y);
    nor = @(x)sqrt(x'*x);
end
if isfield(options,'precd')
    P = options.precd;
else
    P = @(x)x;
end
if isfield(options,'k')
elseif isfield(options,'V')
    options.k = size(options.V,2);
else
    options.k = 1;
    disp('HiOSD: k is reset as 1.');
end
k = options.k;

V = zeros(n,k);
if isfield(options,'V')
    if rank(options.V) < k
        V = randn(n,k);
    else
        V = options.V;
    end
else
    V = randn(n,k);
end

for i = 1 : k
    V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
    V(:,i) = V(:,i)/nor(V(:,i));
end
l=1e-6;
if isfield(options,'minl')
    minl=options.minl;
else
    minl=1e-6;
end
if isfield(options,'maxiter')
    maxiter=options.maxiter;
else
    maxiter=1e3;
end

%% setting environment
alpha (1 : k, 1) = 0;
H = @ (xx,vv,ll) Hvv(xx,vv);
U = zeros(n,k);
W = zeros(n,k);
Y = zeros(n,k);
%% bb iterations
iter = 1;
Vp = zeros(n,k);
Up = zeros(n,k);

% (optional) projector config
if isfield(options, 'proj_cfg') && ~isempty(options.proj_cfg)
    proj_cfg = options.proj_cfg;
else
    proj_cfg = struct();
end
if ~isfield(proj_cfg, 'tol') || isempty(proj_cfg.tol)
    proj_cfg.tol = 1e-6;
end

while iter < maxiter

    % ============================================================
    % 0) Build projector at current x:
    %    Px = projection to normal space of translation orbit at x
    % ============================================================
    proj = lb_translation_projector(x, proj_cfg);  % you already have this file
    Px = proj.PN;                                   % normal projector (operator)

    % Projected Hessian operator: Hproj = Px * H * Px
    % Note: keep same signature H(x,v,l)
    Hproj = @(v) Px( H(x, Px(v), l) );

    % Projected preconditioner: Pproj = Px * P * Px
    Pproj = @(v) Px( P( Px(v) ) );

    % Ensure stored blocks are in normal space (numerical drift)
    V  = Px(V);
    Vp = Px(Vp);

    % Re-orthonormalize V inside normal space
    for i = 1:k
        V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
        V(:,i) = Px(V(:,i));
        V(:,i) = V(:,i) / nor(V(:,i));
    end

    emp = true(3*k, 1);

    % ============================================================
    % 1) Build U, W, Y with projected Hessian and projected precond
    % ============================================================
    for i = 1:k
        % U = Hproj(V)
        U(:,i) = Hproj(V(:,i));

        alpha(i) = innp(V(:,i), U(:,i));

        % Residual direction
        W(:,i) = U(:,i) - alpha(i) * V(:,i);
        W(:,i) = Px(W(:,i));

        % Precondition and re-project
        W(:,i) = Pproj(W(:,i));

        % Orthogonalize against current V and previous W's (all in normal space)
        W(:,i) = W(:,i) - V * innp(V, W(:,i));
        W(:,i) = W(:,i) - W(:,1:i-1) * innp(W(:,1:i-1), W(:,i));
        W(:,i) = Px(W(:,i));

        nrmW = nor(W(:,i));
        if nrmW > 1e-8
            W(:,i) = W(:,i) / nrmW;

            % Y = Hproj(W)
            Y(:,i) = Hproj(W(:,i));
        else
            W(:,i) = 0;
            Y(:,i) = 0;
            emp(k+i) = false;
        end
    end

    % ============================================================
    % 2) Build Vp, Up (previous directions) inside normal space
    % ============================================================
    for i = 1:k
        % Make sure Vp(:,i) stays in normal space before orthogonalization
        Vp(:,i) = Px(Vp(:,i));

        % Orthogonalize against span{V, W, previous Vp}
        B = [V W Vp(:,1:i-1)];
        Vp(:,i) = Vp(:,i) - B * innp(B, Vp(:,i));
        Vp(:,i) = Px(Vp(:,i));

        nrmVpi = nor(Vp(:,i));
        if nrmVpi > 1e-8
            Vp(:,i) = Vp(:,i) / nrmVpi;

            % one more cleanup orthogonalization (like your original)
            B = [V W Vp(:,1:i-1)];
            Vp(:,i) = Vp(:,i) - B * innp(B, Vp(:,i));
            Vp(:,i) = Px(Vp(:,i));
            Vp(:,i) = Vp(:,i) / nor(Vp(:,i));

            % Up = Hproj(Vp)
            Up(:,i) = Hproj(Vp(:,i));
        else
            Vp(:,i) = 0;
            Up(:,i) = 0;
            emp(2*k+i) = false;
        end
    end

    % ============================================================
    % 3) Rayleigh–Ritz on the projected subspace (UU inside normal)
    % ============================================================
    UU = [V, W, Vp];  UU = UU(:, emp);
    YY = [U, Y, Up];  YY = YY(:, emp);

    % Ensure (optional safety): UU,YY are in normal space
    UU = Px(UU);
    YY = Px(YY);

    Pn = innp(UU, YY);
    Pn = (Pn + Pn') / 2;

    [eta, alpha] = eig(Pn, 'vector');
    [alpha, ord] = sort(real(alpha), 'ascend');
    eta = eta(:, ord);

    % Save previous blocks
    Vp = V;
    Up = U;

    % Ritz vectors (then project back to normal, and orthonormalize)
    V = UU * eta(:, 1:k);
    V = Px(V);

    for i = 1:k
        V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
        V(:,i) = Px(V(:,i));
        V(:,i) = V(:,i) / nor(V(:,i));
    end

    alpha = alpha(1:k);

    if mod(iter, 10) == 0
        disp(num2str([iter alpha']));
    end

    iter = iter + 1;
end

end
