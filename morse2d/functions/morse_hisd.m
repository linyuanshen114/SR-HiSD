function [q, V, res, iter] = morse_hisd(q, V, dt , maxiter, a, r0, options)
%MORSE_HISD BB2 position step and LOBPCG direction update. No global state or file I/O.
    tol = options.tol;
    [m, k] = size(V);
    iter = 0;
    
    bb_max_factor = options.bb_max_factor;
    bb_min_factor = options.bb_min_factor;
    Vp = zeros(m,k);
    alph = zeros(k,1);
    W    = zeros(m,k);
    Y    = zeros(m,k);  % Y = HP*W
    [~, g, H] = morse_energy_gradient_hessian(q, a, r0);
    res = norm(g)/sqrt(m);
    
    % Reflect the gradient along the current ascent directions.
    g_proj = V * (V' * g);
    D_curr = g - 2 * g_proj;
    
    q_old = q;
    D_old = D_curr;
    
    q = q - 1e-4 * D_curr; 
    [~, g, H] = morse_energy_gradient_hessian(q, a, r0);
    res = norm(g)/sqrt(m);
    while res > tol && iter < maxiter
        iter = iter + 1;

        g_proj = V * (V' * g); 
        
        D_new = g - 2 * g_proj;
        
        res = norm(g) / sqrt(m);
        
        % BB2 step from the previous displacement and reflected force.
        s_vec = q - q_old;
        y_vec = D_new - D_old;
        sy = s_vec' * y_vec;
        yy = y_vec' * y_vec;
        if yy > eps && isfinite(yy) && isfinite(sy)
            bta = abs(sy / yy);
        else
            bta = dt;
        end
        alpha = min(bta, bb_max_factor * dt);
        alpha = max(alpha, bb_min_factor * dt);
        
        q_old = q;
        D_old = D_new;
        
        q = q - alpha * D_new;

        [~, g, H] = morse_energy_gradient_hessian(q, a, r0);
        res = norm(g)/sqrt(m);

        HP = H;
        [V,~] = qr(V,0);

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
                Y(:,i) = HP * W(:,i);            % Y(:,i) = A W(:,i)
            else
                W(:,i) = 0;
                Y(:,i) = 0;
            end
        end
        
        Up = zeros(m,k);
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
        
        % LOBPCG: Rayleigh--Ritz in direction/residual/previous blocks.
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
        [V,~] = qr(V, 0);
        
    end

end
