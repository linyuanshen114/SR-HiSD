function [q, V, res, iter] = morse3d_srhisd(q, V, dt , maxiter, a, r0, tolerance)
%MORSE3D_SRHISD SR-HiSD in the instantaneous rigid-motion quotient with BB2/LOBPCG.
    [m, k] = size(V);
    res = 1;
    iter = 0;
    alpha = dt;
    Vp = zeros(m-6,k);
    alph = zeros(k,1);
    W    = zeros(m-6,k);
    Y    = zeros(m-6,k);
    Up  = zeros(m-6,k);
    [VMorse, g, H] = morse3d_energy(q, a, r0);
    g_proj = V * (V' * g);
    D_curr = g - 2 * g_proj;
    q_old = q;
    D_old = D_curr;
    % Original warm position step before BB secant updates.
    q = q - 1e-4 * D_curr;
    [VMorse, g, H] = morse3d_energy(q, a, r0);
    while res > tolerance && iter < maxiter
        iter = iter + 1;
        g_proj = V * (V' * g);
        D_new = g - 2 * g_proj;
        res = norm(g) / sqrt(m);
        % BB2 step, clipped to [0.001*dt, 2*dt].
        s_vec = q - q_old;
        y_vec = D_new - D_old;
        sy = s_vec' * y_vec;
        yy = y_vec' * y_vec;
        bta = abs ( sy / yy );
        alpha = min ( bta, 2*dt);
        alpha = max ( alpha, 1e-3*dt);
        q_old = q;
        D_old = D_new;
            q = q - alpha * D_new;
        [VMorse, g, H] = morse3d_energy(q, a, r0);
        % Rebuild the rigid-motion quotient at the new configuration.
        [~, ~,Q] = morse3d_rigid_basis(q);
        HP = Q' * H * Q;
        HP = (HP + HP') / 2;
        V  = Q'*V;
        [V,~] = qr(V,0);
        % Three-block Ritz/LOBPCG update in the quotient coordinates.
        U = HP * V;
        for i = 1:k
            alph(i) = V(:,i)'*U(:,i);
            w       = U(:,i) - alph(i)*V(:,i);
            w = w - V * (V'*w);
            if i > 1
                w = w - W(:,1:i-1) * (W(:, 1:i-1)' * w);
            end
            nrmw = norm(w);
            if nrmw > 1e-8
                W(:,i) = w / nrmw;
                Y(:,i) = HP * W(:,i);
            else
                W(:,i) = 0;
                Y(:,i) = 0;
            end
        end
        Up = zeros(m-6,k);
        for i = 1:k
            vp = Vp(:,i);
            Z = [V, W, Vp(:,1:i-1)];
            vp = vp - Z * (Z'*vp);
            nrmvp = norm(vp);
            if nrmvp > 1e-8
                vp = vp / nrmvp;
                Z  = [V, W, Vp(:,1:i-1)];
                vp = vp - Z * (Z'* vp);
                vp = vp / norm(vp);
                Vp(:,i) = vp;
                Up(:,i) = HP * vp;
            else
                Vp(:,i) = 0;
                Up(:,i) = 0;
            end
        end
        UU = [V, W, Vp];
        YY = [U, Y, Up];
        keep = vecnorm(UU,2,1) > 1e-8;
        UU   = UU(:, keep);
        YY   = YY(:, keep);
        Pn = UU'*YY;
        Pn = (Pn + Pn') / 2;
        [eta, theta] = eig(Pn, 'vector');
        [theta, idx] = sort(theta, 'ascend');
        eta = eta(:, idx);
        Vp  = V;
        V   = UU * eta(:,1:k);
        [V,~] = qr(Q*V, 0);
        eigen = theta(1:k);
    end
    eig_str = sprintf('%.4f ', eigen);
    fprintf('Iter: %4d | V: %10.6f | res: %.1e |step:%.1e \n|eigen:%s\n', ...
    iter, VMorse, res, alpha, eig_str')
end
