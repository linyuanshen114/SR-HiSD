function c=lb_original_candidate(phi,solver,target,output_dir,label,quick)
%LB_ORIGINAL_CANDIDATE Check force and the original mean-zero Hessian index.
% HEX LBSD uses the historical warm/cold raw-HVP gate and a +/-1e-7 zero band.
% This is finite low-spectrum evidence, not a global spectral enclosure.
force=lb_negative_gradient(phi); residual=norm(force)/sqrt(numel(phi));
gate=struct('verified',false,'status','force_not_converged','index',NaN,'signs',[]);
V=[]; lambda=[];
if residual<=1e-7
    try
        [gate,V,lambda]=lb_original_index_gate(phi,solver.V,output_dir,label,quick);
        if gate.verified, gate.signs=gate.starts{1}.signs; else, gate.signs=[]; end
    catch exception
        gate=struct('verified',false,'status','exception','index',NaN, ...
            'signs',[],'error',exception.message);
        fprintf('ORIGINAL_GATE_EXCEPTION %s %s\n',label,exception.message);
    end
end
c=struct('phi',phi,'out',solver,'residual',residual,'force_max',max(abs(force)), ...
    'energy',lb_energy(phi),'index',gate.index,'full_index',gate.index, ...
    'V_full',V,'lambda_full',lambda,'gate',gate,'target',target, ...
    'valid',gate.verified && gate.index==target, ...
    'spectra_finite',all(isfinite(lambda)) && ~isempty(lambda), ...
    'strict_index_certified',false,'independent_index1_verified',gate.verified && gate.index==1);
fprintf('Candidate %s: raw index=%g target=%d RMS=%.6e gate=%s accepted=%d\n', ...
    label,c.index,target,c.residual,gate.status,c.valid);
end
