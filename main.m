%% Calculate all experiments and automatically export Table 1 and Figures 1--3.
% Press Run. These settings are explicit and do not reuse workspace options.
% For a selected task or step, use run_task; see README.md.
review_task='all';
review_subtask='all';
review_options=struct('Mode','compute','Quick',false,'Visible',false);

%% Run the selected scientific scripts in separate workspaces.
review_root=fileparts(mfilename('fullpath'));
review_old_path=path;
review_reports={};
review_figures={};
addpath(review_root,'-begin');
try
    [review_reports,review_figures]=run_task(review_task,review_subtask,review_options);
catch review_error
    path(review_old_path);
    rethrow(review_error);
end
path(review_old_path);

%% A complete run must export all three main figures.
if strcmp(review_task,'all') && strcmp(review_subtask,'all') && ...
        strcmp(review_options.Mode,'compute') && ~review_options.Quick
    review_passed=cellfun(@(report) strcmp(report.status,'PASS'),review_reports);
    assert(all(review_passed),'review:IncompleteCalculation', ...
        'Some calculations are incomplete. See review_reports and the task logs; available results have been saved.');
    for review_number=1:3
        review_name=sprintf('figure%d.png',review_number);
        review_found=cellfun(@(file) endsWith(file,[filesep review_name]) && isfile(file),review_figures);
        assert(any(review_found),'review:MissingFigure', ...
            'The complete calculation did not produce %s. See review_reports and the task logs.',review_name);
    end
end
