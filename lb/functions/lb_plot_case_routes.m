function [outputs,records]=lb_plot_case_routes(routes,output_dir,visible)
%LB_PLOT_CASE_ROUTES Plot actual case states, labeling any missing GF ends.
assert(~isempty(routes),'LBReview:EmptyPlot','Select a case.');
case_id=routes(1).case_id;
assert(all(strcmp({routes.case_id},case_id)),'LBReview:MixedCases','Select one case.');
ready=[routes.complete]; if isfield(routes,'plot_ready'), ready=[routes.plot_ready]; end
assert(any(ready),'LBReview:IncompletePlot','No force-converged state is ready for plotting.');
style=lb_plot_style(case_id); visibility='off'; if visible, visibility='on'; end
outputs={}; records=struct([]); kicks=unique([routes.kick],'stable');
for kick=kicks
    selected=find(abs([routes.kick]-kick)<1e-12); group=routes(selected);
    assert(numel(group)==2 && strcmp(group(2).method,'srhisd1'), ...
        'LBReview:CaseMethods','Expected LBSD and SR-HiSD for each kick.');
    if ~any(ready(selected)), continue; end
    if strcmp(case_id,'bcc2hex')
        stem='transiti1'; caption='BCC to HEX';
    else
        stem='transiti2'; caption='BCC to disorder';
        if abs(kick-.02)>1e-12, stem=sprintf('%s_kick_%04d',stem,round(1000*kick)); end
    end
    columns=max(arrayfun(@(r) numel(r.stages)+2,group));
    fig=figure('Visible',visibility,'Color','w','Position',[60 60 1800 850]);
    annotation(fig,'textbox',[.02 .95 .96 .04], ...
        'String',sprintf('%s | kick %.2g | LBSD (top), SR-HiSD (bottom)',caption,kick), ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',20, ...
        'FontWeight','bold','Color',style.colors.text,'Interpreter','none');
    last_image_ax=[];
    for method_index=1:numel(group)
        r=group(method_index);
        if ~ready(selected(method_index))
            position=lb_panel_position(fig,[.025 .09 .89 .80],columns,1+(method_index-1)*columns);
            annotation(fig,'textbox',[position(1) position(2) .5 position(4)], ...
                'String','Search incomplete; no converged state supplied', ...
                'EdgeColor','none','FontSize',14,'Color',style.colors.text);
            continue;
        end
        measured=[];
        if r.complete
            assert(numel(r.endpoints)==2 && all([r.endpoints.converged]) && ...
                all(strcmp({r.endpoints.parent_phi_sha256},r.parent_phi_sha256)), ...
                'LBReview:NucleusParent','Nucleus endpoints must belong to this saddle.');
            measured=lb_nucleus_measure(r.saddle,r.bcc,r.other,r.L,1400);
            assert(~isempty(measured.ellipse),'LBReview:Nucleus','No critical-nucleus ellipse measured.');
            nucleus_stem=sprintf('%s_%s',case_id,r.method);
            if strcmp(case_id,'bcc2disorder') && abs(kick-.02)>1e-12
                nucleus_stem=sprintf('%s_kick_%04d',nucleus_stem,round(1000*kick));
            end
            %% Only a same-parent pair of actual GF endpoints defines a nucleus.
            single=figure('Visible',visibility,'Color','w','Position',[100 100 740 620]);
            ax=axes(single); imagesc(ax,measured.u,measured.v,measured.phi_slice);
            axis(ax,'image'); axis(ax,'xy'); axis(ax,'off'); colormap(ax,style.cmap); clim(ax,style.clim);
            draw_nucleus(ax,measured.ellipse);
            bar=colorbar(ax); bar.Ticks=style.ticks; bar.Label.String='\phi';
            bar.Color=style.colors.text; bar.FontSize=13;
            title(ax,sprintf('%s | %s | critical nucleus',caption,r.method),'Interpreter','none','Color',style.colors.text);
            single_files={fullfile(output_dir,[nucleus_stem '.png']),fullfile(output_dir,[nucleus_stem '.pdf'])};
            exportgraphics(single,single_files{1},'Resolution',220);
            exportgraphics(single,single_files{2},'ContentType','image','Resolution',300);
            if ~visible, close(single); end
            outputs=[outputs single_files]; %#ok<AGROW>
            record=struct('case_id',case_id,'method',r.method,'kick',kick, ...
                'parent_phi_sha256',r.parent_phi_sha256,'ellipse',measured.ellipse, ...
                'color_limits',style.clim,'colorbar_ticks',style.ticks,'files',{single_files});
            if isempty(records), records=record; else, records(end+1)=record; end %#ok<AGROW>
        end

        %% Plot only supplied fields; an empty GF slot is an explicit text box.
        fields=[{r.mother_bcc},r.stages];
        labels=[{'BCC'},arrayfun(@(k) sprintf('%d-LBSD',k),r.stage_indices,'UniformOutput',false)];
        if method_index==2, labels{2}='1-SR-HiSD'; end
        if r.complete, fields{end+1}=r.other; labels{end+1}=r.other_label; end
        slots=1:numel(fields);
        if columns==6 && method_index==2
            slots=[1 3]; if r.complete, slots=[slots 6]; end
        end
        row_axes=gobjects(1,numel(fields));
        for panel=1:numel(fields)
            ax=axes('Parent',fig,'Units','normalized','Position', ...
                lb_panel_position(fig,[.025 .09 .89 .80],columns,slots(panel)+(method_index-1)*columns));
            row_axes(panel)=ax; is_saddle=panel==numel(r.stages)+1;
            if is_saddle && r.complete
                slice=measured.phi_slice; u=measured.u; v=measured.v;
            else
                [slice,u,v]=lb_oblique_slice(fields{panel},r.L,256);
            end
            imagesc(ax,u,v,slice); axis(ax,'image'); axis(ax,'xy'); axis(ax,'off');
            xlim(ax,[-r.L/2 r.L/2]); ylim(ax,[-r.L/2 r.L/2]);
            colormap(ax,style.cmap); clim(ax,style.clim);
            title(ax,labels{panel},'FontSize',15,'Interpreter','none','Color',style.colors.text);
            if is_saddle && r.complete, draw_nucleus(ax,measured.ellipse); end
        end
        last_image_ax=row_axes(end);
        if ~r.complete
            missing_position=lb_panel_position(fig,[.025 .09 .89 .80],columns,columns+(method_index-1)*columns);
            missing=axes('Parent',fig,'Position',missing_position); axis(missing,'off');
            note=sprintf('GF endpoints\nnot supplied');
            if isfield(r,'connection_note'), note=r.connection_note; end
            text(missing,.5,.5,note,'HorizontalAlignment','center', ...
                'FontSize',16,'Color',style.colors.text,'Tag','missing_gf');
        end
        drawnow;
        for panel=1:numel(row_axes)-1
            left=row_axes(panel).Position; right=row_axes(panel+1).Position;
            arrow_color=style.colors.LBSD; tag='search_arrow';
            if method_index==2, arrow_color=style.colors.SRHiSD; end
            if r.complete && panel==numel(row_axes)-1, arrow_color=style.colors.GF; tag='gf_arrow'; end
            annotation(fig,'arrow',[left(1)+left(3)+.003,right(1)-.003], ...
                [left(2)+left(4)/2,right(2)+right(4)/2],'Color',arrow_color, ...
                'LineWidth',1.3,'HeadLength',6,'HeadWidth',6,'Tag',tag);
        end
    end
    original_position=last_image_ax.Position;
    bar=colorbar(last_image_ax,'Position',[.946 .12 .012 .73]); bar.Ticks=style.ticks;
    bar.Label.String='\phi'; bar.Color=style.colors.text; bar.FontSize=13;
    last_image_ax.Position=original_position;
    assert(numel(findobj(fig,'Tag','critical_nucleus'))==sum([group.complete]));
    assert(numel(findall(fig,'Tag','gf_arrow'))==sum([group.complete]));
    panel_files={fullfile(output_dir,[stem '.png']),fullfile(output_dir,[stem '.pdf'])};
    exportgraphics(fig,panel_files{1},'Resolution',220);
    exportgraphics(fig,panel_files{2},'ContentType','image','Resolution',300);
    if ~visible, close(fig); end
    outputs=[outputs panel_files]; %#ok<AGROW>
end
end

function draw_nucleus(ax,e)
theta=linspace(0,2*pi,400);
rotation=[cos(e.angle) -sin(e.angle);sin(e.angle) cos(e.angle)];
xy=rotation*[e.width/2*cos(theta);e.height/2*sin(theta)];
hold(ax,'on'); plot(ax,e.cx+xy(1,:),e.cy+xy(2,:),'w','LineWidth',2.2,'Tag','critical_nucleus');
end
