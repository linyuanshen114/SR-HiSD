function [V, g, H] = morse3d_energy(q, a, r0)
%MORSE3D_ENERGY Morse energy, analytic gradient and Hessian in 3D.
    if nargin < 3, r0 = 1; end
    if nargin < 2, a = 6;  end
    num_dof = length(q);
    if mod(num_dof,3) ~= 0
        error('q must have 3N coordinates.');
    end
    N = num_dof / 3;
    pos = reshape(q, 3, N)';
    V = 0;
    g = zeros(num_dof,1);
    H = zeros(num_dof,num_dof);
    I3 = eye(3);
    for i = 1:N-1
        idx_i = (i-1)*3+1:i*3;
        pi = pos(i,:)';
        for j = i+1:N
            idx_j = (j-1)*3+1:j*3;
            pj = pos(j,:)';
            dvec = pi - pj;
            r2 = dvec.'*dvec;
            r  = sqrt(r2);
            if r < 1e-12
                r = 1e-12;
                r2 = r^2;
            end
            t  = exp(-a*(r-r0));
            t2 = t*t;
            V = V + (t2 - 2*t);
            dV_dr = 2*a*(t - t2);
            g_scalar = dV_dr / r;
            grad_term = g_scalar * dvec;
            g(idx_i) = g(idx_i) + grad_term;
            g(idx_j) = g(idx_j) - grad_term;
            if nargout > 2
                d2V_dr2 = 2*a*a*(2*t2 - t);
                A = (d2V_dr2 - g_scalar) / r2;
                B = g_scalar;
                H_pair = A*(dvec*dvec.') + B*I3;
                H(idx_i,idx_i) = H(idx_i,idx_i) + H_pair;
                H(idx_j,idx_j) = H(idx_j,idx_j) + H_pair;
                H(idx_i,idx_j) = H(idx_i,idx_j) - H_pair;
                H(idx_j,idx_i) = H(idx_j,idx_i) - H_pair;
            end
        end
    end
end
