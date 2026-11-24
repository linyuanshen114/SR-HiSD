function ctx = review_setup(script_file,subtask,options)
%REVIEW_SETUP Keep paths, run folders and logging out of scientific scripts.
if nargin<3; options=struct(); end
defaults=struct('Mode','compute','Quick',false,'Visible',false, ...
    'RunId','','InputDir','','OutputDir','');
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(options,names{k}); options.(names{k})=defaults.(names{k}); end
end
options.Mode=validatestring(options.Mode,{'compute','reference','check'});
assert(islogical(options.Quick) && isscalar(options.Quick),'review:Quick','Quick must be logical.');
assert(islogical(options.Visible) && isscalar(options.Visible),'review:Visible','Visible must be logical.');
assert(~options.Quick || strcmp(options.Mode,'compute'),'review:Quick','Quick requires compute mode.');
task_dir=fileparts(fileparts(script_file)); root=fileparts(task_dir);
options=review_options_paths(options,task_dir);
[~,task]=fileparts(task_dir);
if isempty(options.RunId)
    stamp=char(datetime('now','Format','yyyyMMdd_HHmmss'));
    token=char(java.util.UUID.randomUUID());
    options.RunId=[options.Mode '_' stamp '_' token(1:8)];
end
assert(~isempty(regexp(char(options.RunId),'^[A-Za-z0-9_-]+$','once')), ...
    'review:RunId','RunId must contain only letters, digits, underscores and hyphens.');
if isempty(options.OutputDir)
    run_root=fullfile(task_dir,'results',options.RunId);
else
    run_root=char(java.io.File(char(options.OutputDir)).getCanonicalPath());
end
run_dir=fullfile(run_root,subtask);
assert(~is_within(run_root,fullfile(task_dir,'functions')) && ...
    ~is_within(run_root,fullfile(task_dir,'scripts')) && ...
    ~is_within(run_root,fullfile(task_dir,'data')) && ...
    ~is_within(task_dir,run_root),'review:Output','Output overlaps task sources or data.');
assert(~is_within(root,run_root),'review:Output','Output cannot be the review root or an ancestor.');
for review_task={'morse2d','morse3d','Onsager','lb'}
    for review_folder={'functions','scripts','data'}
        protected=fullfile(root,review_task{1},review_folder{1});
        assert(~is_within(run_root,protected) && ~is_within(protected,run_root), ...
            'review:Output','Output overlaps a task source/data folder.');
    end
end
protected=fullfile(root,'common');
assert(~is_within(run_root,protected) && ~is_within(protected,run_root), ...
    'review:Output','Output overlaps the shared setup functions.');
old_path=path; old_pwd=pwd; old_rng=rng; old_visibility=get(groot,'defaultFigureVisible');
old_threads=maxNumCompThreads;
had_LB=ismember('LB',who('global'));
global LB
old_LB=LB;
guard=onCleanup(@() restore(old_path,old_pwd,old_rng,old_visibility,old_threads,had_LB,old_LB));
addpath(fullfile(task_dir,'functions'),fullfile(root,'common'),'-begin');
if strcmp(task,'lb'); maxNumCompThreads(4); end
set(groot,'defaultFigureVisible',onoff(options.Visible));
log_file=fullfile(task_dir,'logs',[options.RunId '_' subtask '.log']);
if ~strcmp(options.Mode,'check')
    assert(~isfolder(run_dir) || isempty_content(run_dir),'review:ExistingRun', ...
        'This subtask already has outputs. Choose a new RunId: %s',run_dir);
    if ~isfolder(run_dir); mkdir(run_dir); end
    if ~isfolder(fileparts(log_file)); mkdir(fileparts(log_file)); end
    diary(log_file);
    fprintf('%s / %s | mode=%s | quick=%d\n',task,subtask,options.Mode,options.Quick);
end
ctx=struct('root',root,'task',task,'task_dir',task_dir,'subtask',subtask, ...
    'data_dir',fullfile(task_dir,'data'),'run_root',run_root,'run_dir',run_dir, ...
    'log_file',log_file,'mode',options.Mode,'quick',options.Quick, ...
    'visible',options.Visible,'options',options,'guard',guard);
end

function flag=isempty_content(folder)
d=dir(folder); flag=~any(~ismember({d.name},{'.','..'}));
end

function flag=is_within(child,parent)
child=char(java.io.File(char(child)).getCanonicalPath());
parent=char(java.io.File(char(parent)).getCanonicalPath());
flag=strcmp(child,parent) || startsWith(child,[parent filesep]);
end

function value=onoff(flag)
if flag; value='on'; else; value='off'; end
end

function restore(p,c,r,v,threads,had_LB,old_LB)
diary('off'); path(p); cd(c); rng(r); set(groot,'defaultFigureVisible',v);
maxNumCompThreads(threads);
global LB
if had_LB; LB=old_LB; else; clear global LB; end
end
