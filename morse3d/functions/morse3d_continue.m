function [q_out, res_out, it_out] = morse3d_continue(algo, q_in, dt, maxiter, it_in, res_in, target_index, a, r0, tolerance)
q_out = q_in;
res_out = res_in;
it_out = it_in;
if it_in < maxiter || ~isfinite(res_in)
    return;
end
for ext = 1:2
    [Etmp, gtmp, Htmp] = morse3d_energy(q_out, a, r0);
    [~, ~, Qt] = morse3d_rigid_basis(q_out);
    Hrt = symm(Qt' * Htmp * Qt);
    et = sort(eig(Hrt), 'ascend');
    negt = sum(et < -1e-6);
    if ~(isfinite(Etmp) && isfinite(norm(gtmp)) && isfinite(res_out))
        break;
    end
    close_enough = (norm(gtmp) < 2e-4) || (res_out < 2e-5) || (negt == target_index && norm(gtmp) < 1e-3);
    if ~close_enough
        break;
    end
    ext_iter = maxiter * (2^ext);
    if strcmp(algo, 'SR-HiSD')
        Vc = morse3d_initial_frame(q_out, a, r0, true);
        k_use = min(max(target_index, 1), size(Vc, 2));
        [q_new, ~, res_new, it_new] = morse3d_srhisd(q_out, Vc(:, 1:k_use), dt, ext_iter, a, r0, tolerance);
    else
        Vc = morse3d_initial_frame(q_out, a, r0, false);
        [q_new, ~, res_new, it_new] = morse3d_hisd(q_out, Vc(:, 1:4), dt, ext_iter, a, r0, tolerance);
    end
    q_out = q_new;
    res_out = res_new;
    it_out = it_out + it_new;
    if isfinite(res_out) && res_out < 5e-8
        break;
    end
end
end
function M = symm(M)
M = (M + M') / 2;
end
