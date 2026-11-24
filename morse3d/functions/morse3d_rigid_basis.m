function [U, P_perp, U_perp] = morse3d_rigid_basis(q)
%MORSE3D_RIGID_BASIS Translation/rotation basis and internal-shape projector.
    num_dof = length(q);
    if mod(num_dof,3) ~= 0
        error('q length must be a multiple of 3.');
    end
    N = num_dof / 3;
    pos = reshape(q,3,N)';
    cm  = mean(pos,1);
    pc  = pos - cm;
    x = pc(:,1); y = pc(:,2); z = pc(:,3);
    v_tx = repmat([1;0;0],N,1);
    v_ty = repmat([0;1;0],N,1);
    v_tz = repmat([0;0;1],N,1);
    v_rx = reshape([zeros(N,1), -z, y].',[],1);
    v_ry = reshape([z, zeros(N,1), -x].',[],1);
    v_rz = reshape([-y, x, zeros(N,1)].',[],1);
    A = [v_tx, v_ty, v_tz, v_rx, v_ry, v_rz];
    if nargout > 2
        [Qall,~] = qr(A,0);
        if size(Qall,2) < 6
            error('Degenerate geometry: rigid basis rank < 6.');
        end
        U = Qall(:,1:6);
        [Qfull,~] = qr(A);
        U_perp = Qfull(:,7:end);
    else
        [U,~] = qr(A,0);
        if size(U,2) < 6
            error('Degenerate geometry: rigid basis rank < 6.');
        end
    end
    if nargout > 1
        P_perp = eye(num_dof) - U*U';
    end
end
