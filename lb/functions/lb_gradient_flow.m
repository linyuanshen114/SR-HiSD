function [phi, hist] = lb_gradient_flow(phi0, opts, cfg)
% lb_gradient_flow(phi0, opts, cfg)
% Semi-implicit gradient flow for LB:
%   (I + dt*Lhigh) phi^{n+1} = phi^n - dt*(tau*phi - gamma/2*phi^2 + 1/6*phi^3)
%
% opts:
%   dt, maxiter, tol, outputEvery
%   store_history (true/false), histEvery
%   tag (string)
%
% hist:
%   slices{t}, energies(t)

global LB
N = LB.N;
phi = phi0(:);

if ~isfield(opts,'tag'); opts.tag = 'GF'; end
if ~isfield(opts,'store_history'); opts.store_history = false; end
if ~isfield(opts,'histEvery'); opts.histEvery = 20; end
if ~isfield(opts,'outputEvery'); opts.outputEvery = 200; end
if ~isfield(opts,'tol'); opts.tol = 1e-8; end
use_nesterov = isfield(opts,'nesterov') && opts.nesterov;
if isfield(opts,'nesterov_mu') && ~isempty(opts.nesterov_mu)
    nesterov_mu = opts.nesterov_mu;
else
    nesterov_mu = 0.3;
end
if isfield(opts,'nesterov_warmup') && ~isempty(opts.nesterov_warmup)
    nesterov_warmup = opts.nesterov_warmup;
else
    nesterov_warmup = 20;
end

hist = struct();
hist.slices = {};
hist.energies = [];
hist.tags = {};

mid = floor(N/2)+1;
phi_prev = phi;

for it = 1:opts.maxiter
    if use_nesterov && it > nesterov_warmup
        phi_eval = phi + nesterov_mu * (phi - phi_prev);
        phi_eval = phi_eval - mean(phi_eval);
    else
        phi_eval = phi;
    end

    % residual based on RMS norm of negative gradient
    F = lb_negative_gradient(phi_eval);
    res = sqrt((F'*F)/numel(phi));

    do_draw = ~isfield(cfg,'draw_fig') || cfg.draw_fig;

    if cfg.verbose && mod(it, opts.outputEvery) == 0
        E = lb_energy(phi);
        if use_nesterov
            fprintf('%s-NAG it=%d  res=%.3e  E=%.6e  mu=%.3g\n', opts.tag, it, res, E, nesterov_mu);
        else
            fprintf('%s it=%d  res=%.3e  E=%.6e\n', opts.tag, it, res, E);
        end
        if do_draw
            drawcam(phi, sprintf('%s_%05d', opts.tag, it), cfg);
        end
    end

    if res < opts.tol
        break;
    end

    % Semi-implicit step
    par = LB.par;
    phi3 = reshape(phi_eval, N,N,N);
    phi3 = phi3 - mean(phi3(:));

    % explicit low/nonlinear gradient part
    grad_low = par.tau*phi3 - (par.gamma/2)*phi3.^2 + (1/6)*phi3.^3;

    rhs = phi3 - opts.dt * grad_low;

    rhs_hat = fftn(rhs);
    denom = LB.ones_hat + opts.dt * LB.Lhigh_hat;
    phi_new = real(ifftn(rhs_hat ./ denom));
    phi_new = phi_new - mean(phi_new(:));

    phi_prev = phi;
    phi = phi_new(:);

    if opts.store_history && mod(it, opts.histEvery) == 0
        % store mid-plane slice only (memory-light)
        slice = reshape(phi, N,N,N);
        hist.slices{end+1} = slice(:,:,mid);
        hist.energies(end+1,1) = lb_energy(phi);
        hist.tags{end+1} = sprintf('%s_%05d', opts.tag, it);
    end
end

% final draw
if ~isfield(cfg,'draw_fig') || cfg.draw_fig
    drawcam(phi, sprintf('%s_final', opts.tag), cfg);
end

end
