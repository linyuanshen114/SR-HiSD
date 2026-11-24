function review_result = review_execute(script_file,options)
%REVIEW_EXECUTE Give every actual MATLAB script its own clean workspace.
options=review_options_paths(options,pwd);
p=path; c=pwd; r=rng; v=get(groot,'defaultFigureVisible');
guard=onCleanup(@() restore_session(p,c,r,v)); %#ok<NASGU>
restoredefaultpath;
addpath(fullfile(fileparts(fileparts(fileparts(script_file))),'common'),'-begin');
review_options=options; %#ok<NASGU>
try
    run(script_file);
catch review_exception
    fprintf(2,'%s\n',getReport(review_exception,'extended','hyperlinks','off'));
    rethrow(review_exception);
end
assert(exist('review_result','var')==1,'review:ScriptResult', ...
    'The script did not return review_result: %s',script_file);
end

function restore_session(p,c,r,v)
diary('off'); path(p); cd(c); rng(r); set(groot,'defaultFigureVisible',v);
end
