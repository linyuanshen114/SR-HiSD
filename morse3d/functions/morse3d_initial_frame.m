function V = morse3d_initial_frame(q,a,r0,reduced)
%MORSE3D_INITIAL_FRAME Lowest modes of H or of the internal-shape Hessian Q'*H*Q.
[~,~,H] = morse3d_energy(q,a,r0);
if reduced
    [~,~,Q]=morse3d_rigid_basis(q);
    [Z,d]=eig(Q'*H*Q,'vector');
    [~,order]=sort(d,'ascend'); V=Q*Z(:,order);
else
    [V,d]=eig(H,'vector');
    [~,order]=sort(d,'ascend'); V=V(:,order);
end
end
