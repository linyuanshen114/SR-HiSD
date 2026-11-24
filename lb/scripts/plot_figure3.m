%% Figure 3: actual LBSD and cached SR states for both transitions.
% White nuclei and gray GF arrows require same-parent converged endpoints.
script_dir=fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(fileparts(script_dir)),'common'));
if ~exist('review_options','var'), review_options=struct('Mode','reference'); end
ctx=review_setup(mfilename('fullpath'),'plot_figure3',review_options);
if strcmp(ctx.mode,'check')
    assert(isfile(fullfile(ctx.data_dir,'arrays.h5')) && isfile(fullfile(ctx.data_dir,'store_manifest.json')));
    review_result=review_finish(ctx,struct('status','PASS','numerical_search_performed',false));
else
    [routes,sources]=lb_input_routes(ctx,review_options);
    completed=all([routes.complete]); ready=[routes.complete];
    if isfield(routes,'plot_ready'), ready=[routes.plot_ready]; end
    overlays=struct([]); state_counts=zeros(1,2); plotted_states=struct([]); missing_gf={};
    if any(ready)
        %% Reuse only nucleus measurements with this exact saddle hash.
        records=struct([]); nucleus_file=lb_find_result(ctx,'nuclei',review_options);
        if ~isempty(nucleus_file)
            saved=load(nucleus_file,'records'); records=saved.records;
        end
        hex_style=lb_plot_style('bcc2hex');
        fig=figure('Visible',get(groot,'defaultFigureVisible'),'Color','w','Position',[50 50 1800 1250]);
        annotation(fig,'textbox',[.02 .968 .96 .027],'String','Landau-Brazovskii transitions from BCC', ...
            'EdgeColor','none','HorizontalAlignment','center','FontSize',21, ...
            'FontWeight','bold','Color',hex_style.colors.text,'Interpreter','none');
        case_ids={'bcc2hex','bcc2disorder'};
        for case_index=1:2
            case_id=case_ids{case_index}; style=lb_plot_style(case_id);
            if case_index==1
                bounds=[.025 .53 .89 .35]; columns=6; title_y=.918;
                caption='(a) BCC to HEX | LBSD 4 to 3 to 2 to 1 (top), SR-HiSD (bottom)';
                bar_position=[.946 .53 .012 .35];
            else
                bounds=[.025 .04 .89 .35]; columns=3; title_y=.425;
                caption='(b) BCC to disorder (kick 0.02) | LBSD (top), SR-HiSD (bottom)';
                bar_position=[.946 .04 .012 .35];
            end
            annotation(fig,'textbox',[.025 title_y .89 .03],'String',caption, ...
                'EdgeColor','none','HorizontalAlignment','center','FontSize',17, ...
                'FontWeight','bold','Color',style.colors.text,'Interpreter','none');
            selected=find(strcmp({routes.case_id},case_id));
            assert(numel(selected)==2 && strcmp(routes(selected(2)).method,'srhisd1'));
            last_image_ax=[];
            for method_index=1:2
                route_index=selected(method_index); r=routes(route_index);
                if ~ready(route_index)
                    position=lb_panel_position(fig,bounds,columns,1+(method_index-1)*columns);
                    annotation(fig,'textbox',[position(1) position(2) .5 position(4)], ...
                        'String','Search incomplete; no converged state supplied', ...
                        'EdgeColor','none','FontSize',13,'Color',style.colors.text);
                    continue;
                end
                fields=[{r.mother_bcc},r.stages];
                labels=[{'BCC'},arrayfun(@(k) sprintf('%d-LBSD',k),r.stage_indices,'UniformOutput',false)];
                if method_index==2, labels{2}='1-SR-HiSD'; end
                if r.complete, fields{end+1}=r.other; labels{end+1}=r.other_label; end
                slots=1:numel(fields);
                if case_index==1 && method_index==2
                    slots=[1 3]; if r.complete, slots=[slots 6]; end
                end
                row_axes=gobjects(1,numel(fields));
                for panel=1:numel(fields)
                    ax=axes('Parent',fig,'Units','normalized', ...
                        'Position',lb_panel_position(fig,bounds,columns,slots(panel)+(method_index-1)*columns));
                    row_axes(panel)=ax; is_saddle=panel==numel(r.stages)+1;
                    if is_saddle && r.complete
                        found=[];
                        if ~isempty(records)
                            found=find(strcmp({records.case_id},r.case_id) & strcmp({records.method},r.method));
                        end
                        if isempty(found)
                            assert(numel(r.endpoints)==2 && all([r.endpoints.converged]) && ...
                                all(strcmp({r.endpoints.parent_phi_sha256},r.parent_phi_sha256)));
                            measured=lb_nucleus_measure(r.saddle,r.bcc,r.other,r.L,1400);
                        else
                            assert(isscalar(found) && strcmp(records(found).parent_phi_sha256,r.parent_phi_sha256), ...
                                'LBReview:StaleNucleus','Nucleus cache belongs to another saddle.');
                            measured=records(found).measured;
                        end
                        u=measured.u; v=measured.v; slice=measured.phi_slice;
                    else
                        [slice,u,v]=lb_oblique_slice(fields{panel},r.L,256);
                    end
                    imagesc(ax,u,v,slice); axis(ax,'image'); axis(ax,'xy'); axis(ax,'off');
                    xlim(ax,[-r.L/2 r.L/2]); ylim(ax,[-r.L/2 r.L/2]);
                    colormap(ax,style.cmap); clim(ax,style.clim);
                    title(ax,labels{panel},'FontSize',15,'Interpreter','none','Color',style.colors.text);
                    if is_saddle && r.complete
                        e=measured.ellipse; assert(~isempty(e)); theta=linspace(0,2*pi,400);
                        rotation=[cos(e.angle) -sin(e.angle);sin(e.angle) cos(e.angle)];
                        xy=rotation*[e.width/2*cos(theta);e.height/2*sin(theta)];
                        hold(ax,'on'); plot(ax,e.cx+xy(1,:),e.cy+xy(2,:),'w','LineWidth',2.2,'Tag','critical_nucleus');
                        overlay=struct('case_id',r.case_id,'method',r.method, ...
                            'parent_phi_sha256',r.parent_phi_sha256,'ellipse',e);
                        if isempty(overlays), overlays=overlay; else, overlays(end+1)=overlay; end
                    end
                end
                last_image_ax=row_axes(end);
                if ~r.complete
                    missing_gf{end+1}=[r.case_id '/' r.method]; %#ok<SAGROW>
                    missing_position=lb_panel_position(fig,bounds,columns,columns+(method_index-1)*columns);
                    missing=axes('Parent',fig,'Position',missing_position); axis(missing,'off');
                    note=sprintf('GF endpoints\nnot supplied');
                    if isfield(r,'connection_note'), note=r.connection_note; end
                    text(missing,.5,.5,note,'HorizontalAlignment','center', ...
                        'FontSize',14,'Color',style.colors.text,'Tag','missing_gf');
                end
                drawnow;
                for panel=1:numel(row_axes)-1
                    left=row_axes(panel).Position; right=row_axes(panel+1).Position;
                    arrow_color=style.colors.LBSD; tag='search_arrow';
                    if method_index==2, arrow_color=style.colors.SRHiSD; end
                    if r.complete && panel==numel(row_axes)-1, arrow_color=style.colors.GF; tag='gf_arrow'; end
                    annotation(fig,'arrow',[left(1)+left(3)+.003,right(1)-.003], ...
                        [left(2)+left(4)/2,right(2)+right(4)/2], ...
                        'Color',arrow_color,'LineWidth',1.3,'HeadLength',6,'HeadWidth',6,'Tag',tag);
                end
                state_counts(case_index)=state_counts(case_index)+numel(fields);
                state=struct('case_id',r.case_id,'method',r.method, ...
                    'parent_phi_sha256',r.parent_phi_sha256,'complete',r.complete, ...
                    'stage_indices',r.stage_indices,'has_nucleus',r.complete,'has_gf_arrow',r.complete);
                if isempty(plotted_states), plotted_states=state; else, plotted_states(end+1)=state; end %#ok<SAGROW>
            end
            if ~isempty(last_image_ax)
                original_position=last_image_ax.Position;
                bar=colorbar(last_image_ax,'Position',bar_position); bar.Ticks=style.ticks;
                bar.Label.String='\phi'; bar.Color=style.colors.text; bar.FontSize=13;
                last_image_ax.Position=original_position;
            end
        end
        assert(numel(findobj(fig,'Tag','critical_nucleus'))==sum([routes.complete]));
        assert(numel(findall(fig,'Tag','gf_arrow'))==sum([routes.complete]));
        assert(numel(findall(fig,'Tag','missing_gf'))==numel(missing_gf));
        exportgraphics(fig,fullfile(ctx.run_dir,'figure3.png'),'Resolution',220);
        exportgraphics(fig,fullfile(ctx.run_dir,'figure3.pdf'),'ContentType','image','Resolution',300);
        if ~ctx.visible, close(fig); end
    end
    save(fullfile(ctx.run_dir,'result.mat'),'sources','overlays','state_counts','completed','plotted_states','missing_gf','-v7');
    status='PASS'; if ~completed, status='INCOMPLETE'; end
    if ~any(ready), status='SKIPPED'; end
    review_result=review_finish(ctx,struct('status',status,'completed',completed, ...
        'contains_cases',{{'bcc2hex','bcc2disorder'}},'nucleus_overlays',numel(overlays), ...
        'state_counts',state_counts,'missing_gf',{missing_gf},'plotted_states',plotted_states, ...
        'numerical_search_performed',false, ...
        'reason','Missing same-parent GF endpoints are labeled and have no GF arrow or nucleus.'));
end
