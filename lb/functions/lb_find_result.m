function file=lb_find_result(ctx,subtask,requested_options)
%LB_FIND_RESULT Standalone plots can select the latest compatible calculation.
lookup=ctx.options;
if nargin>=3 && ~isfield(requested_options,'RunId') && isempty(lookup.InputDir) && isempty(lookup.OutputDir)
    lookup=rmfield(lookup,'RunId');
end
file=review_find_result(ctx.task_dir,subtask,lookup);
end
