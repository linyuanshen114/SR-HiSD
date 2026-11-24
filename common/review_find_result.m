function file = review_find_result(task_dir,subtask,options)
%REVIEW_FIND_RESULT Resolve a saved subtask with explicit-input precedence.
% Returns '' when absent. Explicit InputDir never falls back to another run.
file='';
if isfield(options,'InputDir') && ~isempty(options.InputDir)
    base=char(java.io.File(char(options.InputDir)).getCanonicalPath());
    candidates={fullfile(base,subtask,'result.mat'),fullfile(base,'result.mat')};
elseif isfield(options,'OutputDir') && ~isempty(options.OutputDir)
    candidates={fullfile(char(options.OutputDir),subtask,'result.mat')};
elseif isfield(options,'RunId') && ~isempty(options.RunId)
    candidates={fullfile(task_dir,'results',options.RunId,subtask,'result.mat')};
else
    d=dir(fullfile(task_dir,'results','*',subtask,'result.mat'));
    [~,ix]=sort([d.datenum],'descend'); d=d(ix);
    candidates=arrayfun(@(x) fullfile(x.folder,x.name),d,'UniformOutput',false);
end
for k=1:numel(candidates)
    if ~isfile(candidates{k}); continue; end
    metadata=fullfile(fileparts(candidates{k}),'metadata.json');
    if isfile(metadata)
        info=jsondecode(fileread(metadata));
        if ~strcmp(info.subtask,subtask); continue; end
        if isfield(options,'Mode') && ~strcmp(info.mode,options.Mode); continue; end
        if isfield(options,'Quick') && info.quick~=options.Quick; continue; end
    end
    file=candidates{k}; return;
end
end
