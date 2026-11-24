function file=lb_write_stage(output_dir,name,candidate,timing)
%LB_WRITE_STAGE Save one actual search stage without duplicated solver arrays.
phi=candidate.phi; solver_frame=candidate.out.V; %#ok<NASGU>
record=struct('residual',candidate.residual,'energy',candidate.energy, ...
    'index',candidate.index,'valid',candidate.valid,'timing',timing); %#ok<NASGU>
if isfield(candidate,'V_full'); full_frame=candidate.V_full; else; full_frame=[]; end %#ok<NASGU>
if isfield(candidate,'lambda_full'); lambda_full=candidate.lambda_full; else; lambda_full=[]; end %#ok<NASGU>
if isfield(candidate,'lambda'); lambda_reduced=candidate.lambda; else; lambda_reduced=[]; end %#ok<NASGU>
file=fullfile(output_dir,[name '.mat']);
save(file,'phi','solver_frame','full_frame','lambda_full','lambda_reduced','record','-v7.3');
end
