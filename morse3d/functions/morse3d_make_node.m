function node = morse3d_make_node(q, E, gn, morse_index, eigs_red, algo, kind, depth, parent_id, first_seen_case, dt, s0, sg, iter, res, target_ok, stable_ok, seen_sr, seen_hisd)
if nargin < 18, seen_sr = false; end
if nargin < 19, seen_hisd = false; end
node = struct();
node.id = 0;
node.label = '';
node.kind = kind;
node.algo_primary = algo;
node.seen_minimizer = strcmp(algo, 'minimizer');
node.seen_sr = seen_sr;
node.seen_hisd = seen_hisd;
node.seen_upward = strcmp(algo, 'SR-HiSD') && depth > 1;
node.depth = depth;
node.parent_id = parent_id;
node.first_seen_case = first_seen_case;
node.dt = dt;
node.s0 = s0;
node.sign = sg;
node.morse_index = morse_index;
node.E = E;
node.gn = gn;
node.res = res;
node.iter = iter;
node.eig1 = eigs_red(1);
node.eig2 = eigs_red(min(2, numel(eigs_red)));
node.target_ok = target_ok;
node.stable_ok = stable_ok;
node.q = q;
node.kind = kind;
node.target_index = nan;
if morse_index == 0
    node.label = '';
else
    node.label = '';
end
end
