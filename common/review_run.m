function reports = review_run(task,subtask,options)
%REVIEW_RUN Scheduling only; scientific computation is in each task script.
root=fileparts(fileparts(mfilename('fullpath')));
catalog={ ...
    'morse2d',{'original_methods','npss','plot_table1'}; ...
    'morse3d',{'compare_methods','landscape','plot_figure1'}; ...
    'Onsager',{'search','plot_figure2'}; ...
    'lb',{'bcc2hex','bcc2disorder','nuclei','plot_figure3'}};
if nargin<1; task='all'; end
if nargin<2; subtask='all'; end
if nargin<3; options=struct(); end
if ~isfield(options,'Mode'); options.Mode='compute'; end
if ~isfield(options,'Quick'); options.Quick=false; end
if ~isfield(options,'Visible'); options.Visible=false; end
options=review_options_paths(options,pwd);
if strcmpi(task,'list')
    reports=struct('task',catalog(:,1),'subtasks',catalog(:,2)); disp(reports); return;
end
task=validatestring(task,[{'all'},catalog(:,1)']);
assert(strcmp(subtask,'all') || ~strcmp(task,'all'),'review:Subtask', ...
    'Select a task when requesting a single subtask.');
if ~isfield(options,'RunId') || isempty(options.RunId)
    stamp=char(datetime('now','Format','yyyyMMdd_HHmmss'));
    token=char(java.util.UUID.randomUUID()); options.RunId=[stamp '_' token(1:8)];
end
selected=1:size(catalog,1);
if ~strcmp(task,'all'); selected=find(strcmp(catalog(:,1),task)); end
reports={};
for k=selected
    scripts=catalog{k,2};
    if ~strcmp(subtask,'all'); scripts={validatestring(subtask,scripts)}; end
    for j=1:numel(scripts)
        script_file=fullfile(root,catalog{k,1},'scripts',[scripts{j} '.m']);
        fprintf('\n--- %s / %s ---\n',catalog{k,1},scripts{j});
        reports{end+1}=review_execute(script_file,options); %#ok<AGROW>
    end
end
end
