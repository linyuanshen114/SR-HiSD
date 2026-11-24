%% Plot the numerical 3D landscape as states connected by actual graph edges.
% Each inset is the saved 3D particle configuration, not an illustrative image.
if ~exist('review_options','var'); review_options=struct(); end
if ~isfield(review_options,'Mode')
    review_options.Mode='reference';
    if isfield(review_options,'Quick') && review_options.Quick; review_options.Mode='compute'; end
end
task_dir=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(task_dir),'common'));
ctx=review_setup(mfilename('fullpath'),'plot_figure1',review_options);

%% Select the numerical source explicitly.
if strcmp(ctx.mode,'reference') || strcmp(ctx.mode,'check')
    source=fullfile(ctx.data_dir,'reference','result.mat');
    source_kind='saved_validated_result';
else
    lookup=ctx.options;
    if ~isfield(review_options,'RunId') && isfield(lookup,'RunId'); lookup=rmfield(lookup,'RunId'); end
    source=review_find_result(ctx.task_dir,'landscape',lookup);
    assert(~isempty(source),'Morse3D:MissingLandscape','Run landscape.m before drawing a fresh graph.');
    source_kind='computed_landscape';
end
S=load(source,'Nodes','Edges','r0'); Nodes=S.Nodes; Edges=S.Edges;
r0=1; if isfield(S,'r0'); r0=S.r0; end
assert(~isempty(Nodes) && all(isfinite([Nodes.E])));
if strcmp(ctx.mode,'check')
    review_result=review_finish(ctx,struct('status','PASS','nodes',numel(Nodes),'edges',numel(Edges)));
    return;
end

%% Place each state in its Morse-index layer; order within a layer by energy.
position=zeros(numel(Nodes),2); indices=[Nodes.morse_index];
levels=unique(indices); level_count=max(1,numel(levels));
for level=1:numel(levels)
    ids=find(indices==levels(level)); [~,order]=sort([Nodes(ids).E]); ids=ids(order);
    if numel(ids)==1; xs=.5; else; xs=linspace(.13,.87,numel(ids)); end
    y=.16+.65*(level-1)/max(1,level_count-1);
    position(ids,:)=[xs(:),repmat(y,numel(ids),1)];
end
visibility='off'; if ctx.visible; visibility='on'; end
fig=figure('Color','w','Visible',visibility,'Position',[60 60 1650 1050]);
annotation(fig,'textbox',[.05 .94 .9 .045],'String', ...
    sprintf('3D Morse landscape | N = 8 | %d states, %d recorded transitions',numel(Nodes),numel(Edges)), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',20,'FontWeight','bold','Color',[.1 .15 .2]);
colors=struct('SR',[.74 .21 .29],'HiSD',[0 .46 .67],'GD',[.44 .47 .49]);

%% Draw one arrow for each distinct observed source/target/method connection.
seen=cell(0,1); displayed_edges=0;
for j=1:numel(Edges)
    edge=Edges(j); from=edge.source_id; to=edge.target_id;
    key=sprintf('%d/%d/%s',from,to,edge.algo);
    if from==to || any(strcmp(seen,key)); continue; end
    seen{end+1}=key; displayed_edges=displayed_edges+1;
    p1=position(from,:); p2=position(to,:); delta=p2-p1;
    p1=p1+.06*delta/norm(delta); p2=p2-.06*delta/norm(delta);
    if strcmp(edge.algo,'SR-HiSD'); color=colors.SR;
    elseif strcmp(edge.algo,'HiSD'); color=colors.HiSD;
    else; color=colors.GD; end
    annotation(fig,'arrow',[p1(1) p2(1)],[p1(2) p2(2)], ...
        'Color',color,'LineWidth',1,'HeadLength',5,'HeadWidth',5);
end

%% Draw the atoms/bonds and label their physical index and energy.
for j=1:numel(Nodes)
    center=position(j,:);
    ax=axes(fig,'Position',[center(1)-.045 center(2)-.035 .09 .105]);
    morse3d_draw_cluster(ax,Nodes(j).q,r0);
    annotation(fig,'textbox',[center(1)-.058 center(2)-.075 .116 .052], ...
        'String',sprintf('%s | index %d\nE = %.6f',Nodes(j).label,Nodes(j).morse_index,Nodes(j).E), ...
        'EdgeColor','none','BackgroundColor','w','HorizontalAlignment','center','FontSize',9.5,'Color',[.1 .15 .2]);
end
names={'SR-HiSD','HiSD','Gradient descent'}; palette=[colors.SR;colors.HiSD;colors.GD];
for j=1:3
    x=.20+.25*(j-1);
    annotation(fig,'line',[x x+.04],[.055 .055],'Color',palette(j,:),'LineWidth',2);
    annotation(fig,'textbox',[x+.045 .035 .2 .035],'String',names{j}, ...
        'EdgeColor','none','FontSize',11,'Color',[.1 .15 .2]);
end

%% Export the figure and retain the source association in the result metadata.
exportgraphics(fig,fullfile(ctx.run_dir,'figure1.png'),'Resolution',220,'Padding',20);
exportgraphics(fig,fullfile(ctx.run_dir,'figure1.pdf'),'ContentType','image','Resolution',300,'Padding',20);
if ~ctx.visible; close(fig); end

%% Export one labeled particle configuration for every landscape state.
node_images_dir=fullfile(ctx.run_dir,'node_images');
mkdir(node_images_dir);
node_image_files=cell(numel(Nodes),1);
for j=1:numel(Nodes)
    node=Nodes(j);
    node_fig=figure('Color','w','Visible','off','Position',[50 50 600 600], ...
        'PaperUnits','inches','PaperPosition',[0 0 6 6],'PaperSize',[6 6]);
    node_ax=axes(node_fig,'Position',[.04 .04 .92 .82]);
    morse3d_draw_cluster(node_ax,node.q,r0);
    camzoom(node_ax,1.8);
    % The same -135/28 view, blue atoms and gray bonds as the landscape.
    set(findobj(node_ax,'Type','line'),'LineWidth',7);
    set(findobj(node_ax,'Type','surface'),'FaceColor',[.31 .54 .77]);
    annotation(node_fig,'textbox',[.05 .88 .9 .10], ...
        'String',sprintf('%s | index %d\nE = %.6f',node.label,node.morse_index,node.E), ...
        'Interpreter','none','EdgeColor','none','HorizontalAlignment','center', ...
        'FontSize',16,'Color',[.1 .15 .2]);
    node_image_files{j}=fullfile(node_images_dir, ...
        sprintf('node_%02d_%s.png',node.id,node.label));
    print(node_fig,node_image_files{j},'-dpng','-r150');
    close(node_fig);
end
result_file=fullfile(ctx.run_dir,'result.mat');
save(result_file,'source','source_kind','position','displayed_edges', ...
    'node_images_dir','node_image_files','-v7');
review_result=review_finish(ctx,struct('status','PASS','result_file',result_file, ...
    'source',source,'source_kind',source_kind,'nodes',numel(Nodes),'edges',numel(Edges), ...
    'numerical_search_performed',false,'figure',fullfile(ctx.run_dir,'figure1.png'), ...
    'node_images_dir',node_images_dir,'node_image_count',numel(node_image_files)));
