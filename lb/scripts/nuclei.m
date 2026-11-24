%% Critical nuclei: actual same-parent endpoint contrast and PCA ellipses
script_dir=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(script_dir)),'common'));
if ~exist('review_options','var'), review_options=struct('Mode','reference'); end
ctx=review_setup(mfilename('fullpath'),'nuclei',review_options);
if strcmp(ctx.mode,'check')
    assert(isfile(fullfile(ctx.data_dir,'arrays.h5')));
    review_result=review_finish(ctx,struct('status','PASS','numerical_search_performed',false));
else
    [routes,sources]=lb_input_routes(ctx,review_options);
    records=struct([]); completed=all([routes.complete]);
    available=find([routes.complete]);
    if ~isempty(available)
        resolution=1400; if ctx.quick, resolution=128; end
        for j=available
            r=routes(j);
            assert(numel(r.endpoints)==2 && all([r.endpoints.converged]));
            assert(all(strcmp({r.endpoints.parent_phi_sha256},r.parent_phi_sha256)));
            %% rho = K*(TS-BCC)^2 / [K*(TS-BCC)^2 + K*(TS-other)^2]
            % K is the periodic compact mollifier, R=2*sqrt(2)*pi.
            % The helper also performs endpoint alignment and the exact
            % periodic cubic slice interpolation used in the source recipe.
            measured=lb_nucleus_measure(r.saddle,r.bcc,r.other,r.L,resolution);
            % The weighted-PCA ellipse uses rho>0.001 in the central focus
            % disk. Full diameters are 4*sqrt(eigenvalues(covariance)).
            e=measured.ellipse; assert(~isempty(e),'No critical-nucleus ellipse measured.');
            %% Image and white outline use the exact same physical coordinates
            style=lb_plot_style(r.case_id);
            fig=figure('Visible',get(groot,'defaultFigureVisible'),'Color','w','Position',[100 100 740 620]);
            ax=axes(fig); imagesc(ax,measured.u,measured.v,measured.phi_slice);
            axis(ax,'image'); axis(ax,'xy'); axis(ax,'off'); colormap(ax,style.cmap); clim(ax,style.clim);
            hold(ax,'on'); theta=linspace(0,2*pi,400);
            R=[cos(e.angle) -sin(e.angle);sin(e.angle) cos(e.angle)];
            xy=R*[e.width/2*cos(theta); e.height/2*sin(theta)];
            plot(ax,e.cx+xy(1,:),e.cy+xy(2,:),'w','LineWidth',2.2);
            bar=colorbar(ax); bar.Ticks=style.ticks; bar.Label.String='\phi';
            bar.Color=style.colors.text; bar.FontSize=13;
            title(ax,sprintf('%s | %s | critical nucleus',r.case_id,r.method),'Interpreter','none','Color',style.colors.text);
            stem=sprintf('%s_%s',r.case_id,r.method);
            exportgraphics(fig,fullfile(ctx.run_dir,[stem '.png']),'Resolution',220);
            exportgraphics(fig,fullfile(ctx.run_dir,[stem '.pdf']),'ContentType','image','Resolution',300);
            if ~ctx.visible, close(fig); end
            record=struct('case_id',r.case_id,'method',r.method,'kick',r.kick, ...
                'parent_phi_sha256',r.parent_phi_sha256,'measured',measured);
            if isempty(records), records=record; else, records(end+1)=record; end
        end
    end
    save(fullfile(ctx.run_dir,'result.mat'),'records','sources','completed','-v7.3');
    status='PASS';
    if ~completed, status='INCOMPLETE'; end
    if isempty(records), status='SKIPPED'; end
    review_result=review_finish(ctx,struct('status',status,'completed',completed, ...
        'measured_saddles',numel(records),'numerical_search_performed',false, ...
        'missing_endpoint_methods',{{routes(~[routes.complete]).method}}, ...
        'reason','Nuclei require two converged GF endpoints with the same parent; missing pairs are left unmeasured.'));
end
