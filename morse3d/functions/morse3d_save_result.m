function filename = morse3d_save_result(folder,S)
%MORSE3D_SAVE_RESULT Keep MAT/CSV file handling outside the numerical scripts.
filename=fullfile(folder,'result.mat');
save(filename,'-struct','S','-v7');
if isfield(S,'CompareRows') && ~isempty(S.CompareRows)
    flat=rmfield(S.CompareRows,{'node','source_q0','eigs_red','q'});
    writetable(struct2table(flat),fullfile(folder,'compare_rows.csv'));
end
if isfield(S,'Nodes') && ~isempty(S.Nodes)
    writetable(struct2table(rmfield(S.Nodes,'q')),fullfile(folder,'landscape_nodes.csv'));
end
if isfield(S,'Edges') && ~isempty(S.Edges)
    writetable(struct2table(S.Edges),fullfile(folder,'landscape_edges.csv'));
end
end
