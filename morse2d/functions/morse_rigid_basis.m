function [U, P_perp, U_perp] = morse_rigid_basis(q)
%MORSE_RIGID_BASIS Two translations, in-plane rotation, projector and QR complement.
    nd = length(q);
    if mod(nd,2) ~= 0, error('q length must be multiple of 2'); end
    N = nd/2;

    pos = reshape(q,2,N)';
    cm = mean(pos,1);
    pc = pos-cm;
    x = pc(:,1); y = pc(:,2);

    v_tx = repmat([1;0],N,1);
    v_ty = repmat([0;1],N,1);
    v_rz = reshape([-y, x].',[],1);

    A=[v_tx,v_ty,v_rz];

    if nargout>2
        [Qall,~]=qr(A);
        U=Qall(:,1:3);
        U_perp=Qall(:,4:end);
    else
        [U,~]=qr(A,0);
    end

    if nargout>1
        P_perp = eye(nd)-U*U';
    end
end
