function [G,Vwarm,lamwarm]=lb_original_index_gate(phi,warmV,outdir,label,quick)
%LB_ORIGINAL_INDEX_GATE Warm/cold raw mean-zero Hessian low-spectrum check.
% The numerical prefix gate does not establish a global spectral enclosure.
% Normal runs retain the archived k=12, 2500-step, two-hour-per-start limits.
if nargin<5, quick=false; end
runs=cell(1,2); cl=cell(1,2); Vwarm=[]; lamwarm=[];
for j=1:2
    o=struct('k',12,'maxiter',2500,'res_tol',1e-8,'seed',20261117+j, ...
        'time_budget',7200,'check_every',50,'verbose',true,'constraint_basis',[],'warm_V',[]);
    if j==1, o.warm_V=warmV; end
    if quick
        o.k=min(o.k,numel(phi)-2); o.maxiter=4; o.time_budget=30; o.check_every=1;
    end
    fprintf('ORIGINAL_GATE_START %s start=%d\n',label,j);
    sp=lb_raw_spectrum(phi,o); runs{j}=sp;
    save(fullfile(outdir,sprintf('%s_original_start%d.mat',label,j)),'sp','-v7.3');
    cl{j}=lb_classify_zero_band(sp.lambda,sp.raw_residual_rms,1e-7);
    if startsWith(sp.status,'failed:'), cl{j}.verified=false; cl{j}.index=NaN; end
    fprintf('ORIGINAL_GATE_END %s start=%d index=%.0f verified=%d seconds=%.3f\n',label,j,cl{j}.index,cl{j}.verified,sp.seconds);
    if j==1, Vwarm=sp.V; lamwarm=sp.lambda; end
end
agree=cl{1}.verified && cl{2}.verified && cl{1}.index==cl{2}.index;
if agree
    p=min(cl{1}.prefix,cl{2}.prefix);
    uncertain=max(2e-7,10*(runs{1}.raw_residual_rms(1:p)+runs{2}.raw_residual_rms(1:p)));
    agree=all(abs(runs{1}.lambda(1:p)-runs{2}.lambda(1:p))<=uncertain) && isequal(cl{1}.signs(1:p),cl{2}.signs(1:p));
end
G=struct('verified',agree,'index',NaN,'starts',{cl},'statuses',{{runs{1}.status,runs{2}.status}}, ...
    'zero_tol',1e-7,'scope','original mean-zero raw HVP; cold/warm numerical low spectrum, not global enclosure');
if agree, G.index=cl{1}.index; G.status='numerically_resolved'; else, G.status='unresolved'; end
save(fullfile(outdir,[label '_original_gate.mat']),'G','-v7.3');
end
