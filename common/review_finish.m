function result = review_finish(ctx,summary)
%REVIEW_FINISH Record a short subtask report and close its log/session guard.
if nargin<2; summary=struct(); end
status='PASS';
if isfield(summary,'status'); status=summary.status; end
figure_files={};
if ~strcmp(ctx.mode,'check')
    files=dir(fullfile(ctx.run_dir,'**','*'));
    files=files(~[files.isdir]);
    for k=1:numel(files)
        [~,~,extension]=fileparts(files(k).name);
        if any(strcmpi(extension,{'.png','.pdf','.svg','.fig','.jpg','.jpeg','.eps'}))
            figure_files{end+1}=fullfile(files(k).folder,files(k).name); %#ok<AGROW>
        end
    end
    figure_files=sort(figure_files);
end
result=struct('task',ctx.task,'subtask',ctx.subtask,'mode',ctx.mode, ...
    'quick',ctx.quick,'status',status,'output_dir',ctx.run_dir, ...
    'result_file',fullfile(ctx.run_dir,'result.mat'),'log_file',ctx.log_file, ...
    'run_id',ctx.options.RunId,'figure_files',{figure_files}, ...
    'outputs',{figure_files},'summary',summary);
if ~strcmp(ctx.mode,'check')
    fid=fopen(fullfile(ctx.run_dir,'metadata.json'),'w'); assert(fid>=0);
    fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true)); fclose(fid);
    fprintf('Finished %s / %s: %s\nResults: %s\n',ctx.task,ctx.subtask,status,ctx.run_dir);
    if ~isempty(figure_files)
        fprintf('Figures (%d):\n',numel(figure_files));
        fprintf('  %s\n',figure_files{:});
    elseif strcmp(status,'SKIPPED') && isfield(summary,'reason')
        fprintf('Plot skipped: %s\n',summary.reason);
    elseif strcmp(status,'INCOMPLETE')
        fprintf('Figures not exported: the numerical search is incomplete.\n');
    end
end
delete(ctx.guard);
end
