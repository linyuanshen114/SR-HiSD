function [q, res, iter] = morse3d_gradient_descent(q, dt, maxiter, a, r0, tolerance)
%MORSE3D_GRADIENT_DESCENT Gradient descent retaining the original returned-residual stopping.
    m = numel(q);
    res = 1;
    iter = 0;
    while res > tolerance && iter < maxiter
        iter = iter + 1;
        [VMorse, g] = morse3d_energy(q, a, r0);
        q = q - dt * g;
        res = norm(g) / sqrt(m);
        if mod(iter, 100) == 0
            fprintf('GD Iter: %4d | V: %10.6f | res: %.1e\n', iter, VMorse, res);
        end
    end
    fprintf('GD Iter: %4d | V: %10.6f | res: %.1e\n', iter, VMorse, res);
end
