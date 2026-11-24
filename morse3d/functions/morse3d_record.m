function out = morse3d_record(q, q0, qmin, a, r0, res, iter, algo, dt, s0, sg, depth, target_index)
[E, g, H] = morse3d_energy(q, a, r0);
[~, ~, Q] = morse3d_rigid_basis(q);
Hred = symm(Q' * H * Q);
eigs_red = sort(eig(Hred), 'ascend');
neg = sum(eigs_red < -1e-6);
stable_ok = isfinite(E) && isfinite(norm(g)) && isfinite(res) && ~any(isnan(q(:))) && ~any(isinf(q(:)));
target_ok = stable_ok && (neg == target_index) && (norm(g) < 1e-5) && (res < 5e-8);
kind = 'saddle';
if neg == 0
    kind = 'minimum';
end
out = struct();
out.algo = algo;
out.dt = dt;
out.s0 = s0;
out.sign = sg;
out.depth = depth;
out.parent_id = 0;
out.target_node_id = 0;
out.target_is_new = false;
out.source_id = 0;
out.E = E;
out.gn = norm(g);
out.res = res;
out.iter = iter;
out.morse_index = neg;
out.eig1 = eigs_red(1);
out.eig2 = eigs_red(min(2, numel(eigs_red)));
out.eigs_red = eigs_red;
out.dist_ref = norm(q - q0);
out.dist_qmin = norm(q - qmin);
out.target_ok = target_ok;
out.target_index = target_index;
out.stable_ok = stable_ok;
out.kind = kind;
out.q = q;
out.node = morse3d_make_node(q,E,norm(g),neg,eigs_red,algo,kind,depth,0,0,dt,s0,sg,iter,res,target_ok,stable_ok,strcmp(algo,'SR-HiSD'),strcmp(algo,'HiSD'));
out.node.target_index = target_index;
out.source_algo = algo;
out.source_q0 = q0;
out.time_sec = NaN;
out.k_target = 0;
end
function M = symm(M)
M = (M + M') / 2;
end
