function [Vnew, Vprevious, info] = lb_spectrum_lobpcg_step(Hv, V, Vprevious, project, precondition, opts)
%LB_SPECTRUM_LOBPCG_STEP Gram-whitened matrix-free low-spectrum Ritz step.
% Inner products are X'*Y/m. The trial basis contains the current vectors,
% projected preconditioned residuals and prior directions. Dependent columns
% are dropped; the small Gram matrix is whitened before the Ritz solve.
% Hv and precondition act on one vector; project also accepts blocks.
% No grid-square Hessian or projector is formed.

if nargin < 3 || isempty(Vprevious), Vprevious = []; end
if nargin < 4 || isempty(project), project = @(X) X; end
if nargin < 5, precondition = []; end
if nargin < 6 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'block_tol'), opts.block_tol = 1e-8; end
if ~isfield(opts, 'gram_tol'), opts.gram_tol = max(100 * eps, opts.block_tol^2); end
assert(isa(Hv, 'function_handle') && isa(project, 'function_handle'), ...
    'lb_spectrum_lobpcg_step:Callbacks', 'Hv and project must be function handles.');
assert(isempty(precondition) || isa(precondition, 'function_handle'), ...
    'lb_spectrum_lobpcg_step:Callbacks', 'precondition must be empty or a function handle.');
assert(isnumeric(V) && ismatrix(V) && ~isempty(V) && all(isfinite(V(:))), ...
    'lb_spectrum_lobpcg_step:InitialBlock', 'V must be a finite, nonempty numeric matrix.');
assert(isscalar(opts.block_tol) && isfinite(opts.block_tol) && opts.block_tol >= 0 ...
    && isscalar(opts.gram_tol) && isfinite(opts.gram_tol) && opts.gram_tol > 0, ...
    'lb_spectrum_lobpcg_step:Options', 'Invalid block_tol or gram_tol.');
[m, k] = size(V);
if isempty(Vprevious), Vprevious = zeros(m, 0, 'like', V); end
assert(isnumeric(Vprevious) && size(Vprevious, 1) == m ...
    && all(isfinite(Vprevious(:))), 'lb_spectrum_lobpcg_step:HistoryBlock', ...
    'Vprevious must have m rows and finite entries, or be empty.');

hvp_count = 0;
precondition_count = 0;
project_count = 0;
[V, initial_kept] = orthogonalize(V, zeros(m, 0, 'like', V));
assert(numel(initial_kept) == k, 'lb_spectrum_lobpcg_step:InitialRank', ...
    'Projected initial V has rank %d, but %d directions are required.', numel(initial_kept), k);
HV = hessian_action(V);

% Rotate within the old block: this leaves its span unchanged and gives
% scalar Ritz values for the individual residual/preconditioner columns.
[rotation, old_values] = eig(hermitian(V' * HV / m), 'vector');
[old_values, order] = sort(real(old_values), 'ascend');
rotation = rotation(:, order);
V = V * rotation;
HV = HV * rotation;
Vold = V;
R = apply_project(HV - V * diag(old_values));
old_residual_rms = sqrt(sum(abs(R).^2, 1) / m);

Wraw = zeros(size(R), 'like', R);
for j = 1:k
    % Avoid calling a user preconditioner on a numerically zero residual.
    % The same absolute RMS threshold is used for candidate rank dropping.
    if old_residual_rms(j) <= opts.block_tol, continue; end
    if isempty(precondition)
        Wraw(:, j) = R(:, j);
    else
        z = precondition(R(:, j));
        precondition_count = precondition_count + 1;
        assert(isequal(size(z), [m, 1]) && all(isfinite(z(:))), ...
            'lb_spectrum_lobpcg_step:PreconditionerOutput', ...
            'precondition must return one finite m-by-1 vector.');
        Wraw(:, j) = z;
    end
end
Wraw = apply_project(Wraw);
[W, ~] = orthogonalize(Wraw, V);
[Pold, ~] = orthogonalize(Vprevious, [V, W]);
B = [V, W, Pold];
HB = [HV, hessian_action(W), hessian_action(Pold)];

G = hermitian(B' * B / m);
T = hermitian(B' * HB / m);
[whitener, gram_rank] = gram_whitener(G);
assert(gram_rank >= k, 'lb_spectrum_lobpcg_step:TrialRank', ...
    'The projected trial space has fewer than k numerically independent directions.');
[small_vectors, values] = eig(hermitian(whitener' * T * whitener), 'vector');
[~, order] = sort(real(values), 'ascend');
coefficients = whitener * small_vectors(:, order(1:k));
Vnew = apply_project(B * coefficients);
HVnew = apply_project(HB * coefficients);

% Small selected-space normalization and Ritz cleanup. Transform the cached
% Hessian images by the SAME coefficients; no extra matrix-free action.
[normalizer, selected_rank] = gram_whitener(hermitian(Vnew' * Vnew / m));
assert(selected_rank == k, 'lb_spectrum_lobpcg_step:SelectedRank', ...
    'The selected Ritz block lost rank.');
Vnew = Vnew * normalizer;
HVnew = HVnew * normalizer;
[rotation, values] = eig(hermitian(Vnew' * HVnew / m), 'vector');
[values, order] = sort(real(values), 'ascend');
rotation = rotation(:, order);
Vnew = Vnew * rotation;
HVnew = HVnew * rotation;
Vprevious = Vold;

residual = HVnew - Vnew * diag(values);
leak = apply_project(Vnew) - Vnew;
info = struct('ritz_values', values(:), ...
    'residual_rms', sqrt(sum(abs(residual).^2, 1) / m).', ...
    'orthogonality_error', norm(Vnew' * Vnew / m - eye(k), 'fro'), ...
    'constraint_leak_rms', sqrt(sum(abs(leak).^2, 1) / m).', ...
    'old_ritz_values', old_values(:), 'old_residual_rms', old_residual_rms(:), ...
    'block_sizes', [k, size(W, 2), size(Pold, 2)], ...
    'trial_dimension', size(B, 2), 'gram_rank', gram_rank, ...
    'hvp_count', hvp_count, 'precondition_count', precondition_count, ...
    'project_count', project_count, 'status', 'ok');

    function Y = apply_project(X)
        if isempty(X), Y = X; return; end
        Y = project(X);
        project_count = project_count + 1;
        assert(isequal(size(Y), size(X)) && all(isfinite(Y(:))), ...
            'lb_spectrum_lobpcg_step:ProjectorOutput', ...
            'project must preserve the block size and return finite values.');
    end

    function Y = hessian_action(X)
        X = apply_project(X);
        Y = zeros(size(X), 'like', X);
        for jj = 1:size(X, 2)
            z = Hv(X(:, jj));
            hvp_count = hvp_count + 1;
            assert(isequal(size(z), [m, 1]) && all(isfinite(z(:))), ...
                'lb_spectrum_lobpcg_step:HessianOutput', ...
                'Hv must return one finite m-by-1 vector.');
            Y(:, jj) = z;
        end
        Y = apply_project(Y);
    end

    function [Q, kept] = orthogonalize(C, against)
        C = apply_project(C);
        Q = zeros(m, 0, 'like', V);
        kept = zeros(1, 0);
        for jj = 1:size(C, 2)
            w = C(:, jj);
            % Two modified Gram-Schmidt passes, with RMS inner products.
            for pass = 1:2
                for ii = 1:size(against, 2)
                    w = w - against(:, ii) * (against(:, ii)' * w / m);
                end
                for ii = 1:size(Q, 2)
                    w = w - Q(:, ii) * (Q(:, ii)' * w / m);
                end
            end
            w = apply_project(w);
            scale = norm(w) / sqrt(m);
            if scale > opts.block_tol
                Q(:, end + 1) = w / scale; %#ok<AGROW>
                kept(end + 1) = jj; %#ok<AGROW>
            end
        end
    end

    function [Z, numerical_rank] = gram_whitener(small_G)
        [U, d] = eig(hermitian(small_G), 'vector');
        d = real(d);
        cutoff = opts.gram_tol * max(1, max(abs(d)));
        retain = d > cutoff;
        numerical_rank = nnz(retain);
        Z = U(:, retain) ./ sqrt(d(retain)).';
    end
end

function A = hermitian(A)
A = (A + A') / 2;
end
