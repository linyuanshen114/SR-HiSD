function [Nodes,id,is_new] = morse3d_merge(Nodes,candidate,unique_tol,minimum_isometry,r0)
if nargin<4; minimum_isometry=false; end
if nargin<5; r0=1; end
if minimum_isometry
    [Nodes,id,is_new]=add_or_merge_minimum_iso(Nodes,candidate,unique_tol,1e-4,2e-2,r0);
else
    [Nodes,id,is_new]=add_or_merge_node(Nodes,candidate,unique_tol);
end
end
function [Nodes, node_id, is_new] = add_or_merge_node(Nodes, candidate, unique_tol)
is_new = true;
node_id = numel(Nodes) + 1;
for ii = 1:numel(Nodes)
    if norm(Nodes(ii).q - candidate.q) < unique_tol && Nodes(ii).morse_index == candidate.morse_index
        node_id = ii;
        is_new = false;
        Nodes(ii) = merge_node(Nodes(ii), candidate);
        return;
    end
end
candidate.id = node_id;
candidate.label = default_label(candidate);
if isempty(Nodes)
    Nodes = candidate;
else
    Nodes(node_id) = candidate; %#ok<AGROW>
end
end
function [Nodes, node_id, is_new] = add_or_merge_minimum_iso(Nodes, candidate, unique_tol, energy_tol, sig_tol, r0)
if strcmp(candidate.kind, 'minimum')
    sig_c = pairwise_signature(candidate.q);
    for ii = 1:numel(Nodes)
        if ~strcmp(Nodes(ii).kind, 'minimum')
            continue;
        end
        if abs(Nodes(ii).E - candidate.E) > energy_tol
            continue;
        end
        sig_i = pairwise_signature(Nodes(ii).q);
        if signatures_close(sig_i, sig_c, sig_tol, r0)
            node_id = ii;
            is_new = false;
            Nodes(ii) = merge_node(Nodes(ii), candidate);
            return;
        end
    end
end
[Nodes, node_id, is_new] = add_or_merge_node(Nodes, candidate, unique_tol);
end
function sig = pairwise_signature(q)
pos = reshape(q, 3, []).';
n = size(pos, 1);
d = zeros(n * (n - 1) / 2, 1);
k = 0;
for i = 1:n-1
    for j = i+1:n
        k = k + 1;
        d(k) = norm(pos(i,:) - pos(j,:));
    end
end
sig = sort(d, 'ascend');
end
function ok = signatures_close(sig1, sig2, tol, r0)
if numel(sig1) ~= numel(sig2)
    ok = false;
    return;
end
scale = max(r0, 1);
err = norm(sig1 - sig2) / max(1, numel(sig1));
ok = err < tol * scale;
end
function node = merge_node(node, candidate)
node.seen_minimizer = node.seen_minimizer || candidate.seen_minimizer;
node.seen_sr = node.seen_sr || candidate.seen_sr;
node.seen_hisd = node.seen_hisd || candidate.seen_hisd;
node.seen_upward = node.seen_upward || candidate.seen_upward;
if strcmp(node.algo_primary, candidate.algo_primary)
    return;
end
if strcmp(node.algo_primary, 'both')
    return;
end
if strcmp(candidate.algo_primary, 'both')
    node.algo_primary = 'both';
elseif strcmp(node.algo_primary, 'minimizer') || strcmp(candidate.algo_primary, 'minimizer')
    node.algo_primary = 'both';
elseif ~strcmp(node.algo_primary, candidate.algo_primary)
    node.algo_primary = 'both';
end
end
function label = default_label(node)
if strcmp(node.kind, 'minimum')
    label = sprintf('m%02d', node.id);
else
    label = sprintf('s%02d', node.id);
end
end
