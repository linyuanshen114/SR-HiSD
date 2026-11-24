function [reports,figures]=run_task(task,subtask,options)
%RUN_TASK Calculate and plot a task, or run a selected step.
%   run_task('morse3d')                  % calculation and Figure 1
%   run_task('Onsager')                  % calculation and Figure 2
%   run_task('lb')                       % both cases, nuclei and Figure 3
%   run_task('all','all',struct('Mode','reference'))
%   run_task('Onsager','search',struct('Mode','compute','RunId','my_run'))
% Outputs stay in <task>/results/<RunId>/<subtask>/ by default.
if nargin<1; task='all'; end
if nargin<2; subtask='all'; end
if nargin<3; options=struct(); end
root=fileparts(mfilename('fullpath'));
old_path=path;
guard=onCleanup(@() path(old_path)); %#ok<NASGU>
addpath(fullfile(root,'common'),'-begin');
reports=review_run(task,subtask,options);
figures={};
if iscell(reports)
    for step=1:numel(reports)
        figures=[figures,reports{step}.figure_files]; %#ok<AGROW>
    end
end
end
