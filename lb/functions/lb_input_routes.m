function [routes,sources]=lb_input_routes(ctx,requested_options)
%LB_INPUT_ROUTES Select the current saved cases or actual same-run results.
if strcmp(ctx.mode,'reference')
    hex=lb_reference_routes(ctx.data_dir,'bcc2hex');
    dis=lb_reference_routes(ctx.data_dir,'bcc2disorder');
    routes=[hex dis];
    sources={fullfile(ctx.data_dir,'arrays.h5'), ...
        fullfile(ctx.data_dir,'reference','bcc2hex_srhisd','final.mat'), ...
        fullfile(ctx.data_dir,'reference','bcc2disorder_srhisd','final.mat')};
    hex_base=fullfile(ctx.data_dir,'reference','bcc2hex_lbsd');
    sources=[sources,{fullfile(ctx.data_dir,'input','bcc2hex_lbsd.mat'), ...
        fullfile(ctx.data_dir,'input','bcc2hex_lbsd_reference.mat'), ...
        fullfile(hex_base,'source_manifest.json'),fullfile(hex_base,'search','status.json'), ...
        fullfile(hex_base,'search','attempts.csv'),fullfile(hex_base,'gf','source.mat'), ...
        fullfile(hex_base,'gf','source.json'),fullfile(hex_base,'gf','status.json')}];
    for label={'up_4','down_4to3_kick_1000_sign_+1', ...
            'down_3to2_kick_1000_sign_+1','down_2to1_kick_1000_sign_-1'}
        sources=[sources,{fullfile(hex_base,'search',[label{1} '_candidate.mat']), ...
            fullfile(hex_base,'search',[label{1} '_original_gate.mat'])}]; %#ok<AGROW>
    end
    for branch={'sign_+1','sign_-1'}
        sources=[sources,{fullfile(hex_base,'gf',branch{1},'initial.mat'), ...
            fullfile(hex_base,'gf',branch{1},'endpoint.mat')}]; %#ok<AGROW>
    end
    for name={'bcc2hex_srhisd','bcc2disorder_srhisd'}
        gf_dir=fullfile(ctx.data_dir,'reference',name{1},'gf');
        if isfile(fullfile(gf_dir,'completion.json'))
            sources=[sources,{fullfile(gf_dir,'completion.json'), ...
                fullfile(gf_dir,'plus','endpoint.mat'),fullfile(gf_dir,'minus','endpoint.mat')}]; %#ok<AGROW>
        end
    end
else
    if nargin<2, requested_options=ctx.options; end
    sources={lb_find_result(ctx,'bcc2hex',requested_options), ...
        lb_find_result(ctx,'bcc2disorder',requested_options)};
    assert(all(cellfun(@isfile,sources)),'LBReview:MissingCalculation', ...
        'Run bcc2hex and bcc2disorder first, or choose reference mode explicitly.');
    parent1=fileparts(fileparts(sources{1})); parent2=fileparts(fileparts(sources{2}));
    assert(strcmp(parent1,parent2),'LBReview:MixedRuns','Both case results must belong to the same run.');
    H=load(sources{1},'routes'); D=load(sources{2},'routes');
    dis=D.routes(abs([D.routes.kick]-.02)<1e-12);
    assert(numel(H.routes)==2 && numel(dis)==2,'LBReview:Routes','Expected both methods for both cases.');
    routes=[H.routes dis];
end
end
