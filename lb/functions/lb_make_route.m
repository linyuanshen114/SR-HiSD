function r=lb_make_route(case_id,method,seed,selected,stages,endpoints,quick,stage_indices)
%LB_MAKE_ROUTE Store numerical fields and explicit completeness/provenance.
bcc=seed.phi_bcc; other=[]; other_label='unconverged';
connection_note=sprintf('GF endpoints\nnot supplied');
negative_verified=isfield(selected,'negative_mode_verified') && selected.negative_mode_verified;
index1_verified=isfield(selected,'independent_index1_verified') && selected.independent_index1_verified;
force_converged=all(isfinite(selected.phi)) && selected.residual<=1.05*seed.bb.epsf;
complete=selected.valid && (selected.index==1 || negative_verified) && ...
    numel(endpoints)==2 && all([endpoints.converged]);
if complete
    is_bcc=strcmpi({endpoints.class},'BCC');
    complete=nnz(is_bcc)==1;
    if complete
        bcc=endpoints(is_bcc).phi; other=endpoints(~is_bcc).phi; other_label=endpoints(~is_bcc).class;
        if strcmp(case_id,'bcc2disorder')
            complete=strcmpi(other_label,'disorder');
        elseif strcmp(case_id,'bcc2hex')
            complete=strcmpi(other_label,'HEX');
        end
        if ~complete, connection_note=sprintf('Endpoint phase\nunresolved'); end
    else
        connection_note=sprintf('BCC connection\nunresolved');
    end
elseif ~isempty(endpoints)
    connection_note=sprintf('GF endpoints\nnot converged');
end
if complete, connection_note=''; end
r=struct('case_id',case_id,'method',method,'N',seed.N,'L',seed.L,'par',seed.par,'kick',seed.kick, ...
    'bcc',bcc,'mother_bcc',seed.phi_bcc,'saddle',selected.phi,'stages',{stages},'stage_indices',stage_indices, ...
    'other',other,'other_label',other_label,'endpoints',endpoints, ...
    'complete',complete,'fresh',true,'quick',quick,'parent_phi_sha256',lb_numeric_sha(selected.phi), ...
    'force_converged',force_converged,'negative_mode_verified',negative_verified, ...
    'strict_index_certified',false,'plot_ready',~quick && selected.valid && force_converged, ...
    'independent_index1_verified',index1_verified,'connection_note',connection_note);
end
