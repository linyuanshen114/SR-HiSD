function report=morse3d_check_paper_result(Nodes,Edges,CompareRows)
%MORSE3D_CHECK_PAPER_RESULT Article-level counts/labels/energies and self-checks.
% No numerical reference field, trajectory or iteration count is compared.
task_dir=fileparts(fileparts(mfilename('fullpath')));
paper=readtable(fullfile(task_dir,'data','paper','figure1_states.csv'),'TextType','string');
report=struct('passed',false,'expected_states',17,'expected_edges',106, ...
    'expected_trials',106,'state_count',numel(Nodes),'edge_count',numel(Edges), ...
    'trial_count',numel(CompareRows),'paper_energy_tolerance',5.1e-7, ...
    'maximum_displayed_energy_error',Inf,'maximum_energy_recompute_error',0, ...
    'maximum_gradient_norm',0,'physical_indices_valid',true, ...
    'trial_certificates_valid',true,'edge_records_valid',true,'connected_states',true);
if numel(Nodes)~=17 || numel(Edges)~=106 || numel(CompareRows)~=106; return; end
assert(height(paper)==17,'Morse3D:PaperData','Figure 1 must list seventeen states.');
if ~isequal([Nodes.id],1:17) || ...
        ~isequal(sort(string({Nodes.label})'),sort(paper.label)) || ...
        ~isequal([Edges.id],1:106)
    return;
end
paper_energy_error=0;
for j=1:numel(Nodes)
    node=Nodes(j); k=find(paper.label==string(node.label));
    if ~isscalar(k) || paper.morse_index(k)~=node.morse_index; return; end
    paper_energy_error=max(paper_energy_error,abs(node.E-paper.energy(k)));
    [energy,g,H]=morse3d_energy(node.q,6,1);
    [~,~,Q]=morse3d_rigid_basis(node.q);
    A=Q'*H*Q; physical_index=sum(eig((A+A')/2)<-1e-6);
    report.maximum_energy_recompute_error=max(report.maximum_energy_recompute_error,abs(energy-node.E));
    report.maximum_gradient_norm=max(report.maximum_gradient_norm,norm(g));
    report.physical_indices_valid=report.physical_indices_valid && physical_index==node.morse_index;
    positions=reshape(node.q,3,[])';
    distance=sqrt(sum((permute(positions,[1 3 2])-permute(positions,[3 1 2])).^2,3));
    report.connected_states=report.connected_states && max(conncomp(graph(distance<=1.6,'omitselfloops')))==1;
end
report.maximum_displayed_energy_error=paper_energy_error;
for j=1:numel(CompareRows)
    row=CompareRows(j); [energy,g,H]=morse3d_energy(row.q,6,1);
    [~,~,Q]=morse3d_rigid_basis(row.q); A=Q'*H*Q;
    physical_index=sum(eig((A+A')/2)<-1e-6);
    accepted=all(isfinite(row.q)) && isfinite(energy) && ...
        physical_index==row.target_index && norm(g)<1e-5 && row.res<5e-8;
    report.trial_certificates_valid=report.trial_certificates_valid && ...
        physical_index==row.morse_index && abs(energy-row.E)<1e-9 && accepted==logical(row.target_ok);
end
for j=1:numel(Edges)
    edge=Edges(j); valid=ismember(edge.source_id,1:17) && ismember(edge.target_id,1:17);
    if valid
        valid=edge.source_index==Nodes(edge.source_id).morse_index && ...
            edge.target_index==Nodes(edge.target_id).morse_index && ...
            abs(edge.source_E-Nodes(edge.source_id).E)<1e-9 && ...
            abs(edge.target_E-Nodes(edge.target_id).E)<1e-9;
    end
    report.edge_records_valid=report.edge_records_valid && valid;
end
report.passed=paper_energy_error<=report.paper_energy_tolerance && ...
    report.maximum_energy_recompute_error<1e-9 && report.maximum_gradient_norm<1e-5 && ...
    report.physical_indices_valid && report.trial_certificates_valid && ...
    report.edge_records_valid && report.connected_states;
end
