function [Edges, edge_id] = morse3d_add_edge(Edges, Nodes, source_id, target_id, row, new_node, depth)
edge_id = numel(Edges) + 1;
edge = struct();
edge.id = edge_id;
edge.source_id = source_id;
edge.target_id = target_id;
edge.algo = row.algo;
edge.depth = depth;
edge.k_target = row.k_target;
edge.dt = row.dt;
edge.s0 = row.s0;
edge.sign = row.sign;
edge.new_node = new_node;
edge.source_E = Nodes(source_id).E;
edge.target_E = Nodes(target_id).E;
edge.source_index = Nodes(source_id).morse_index;
edge.target_index = Nodes(target_id).morse_index;
edge.deltaE = edge.target_E - edge.source_E;
edge.source_label = Nodes(source_id).label;
edge.target_label = Nodes(target_id).label;
if isempty(Edges)
    Edges = edge;
else
    Edges(edge_id) = edge; %#ok<AGROW>
end
end
