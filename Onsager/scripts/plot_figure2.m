%% Draw the four Onsager states. All axes and annotations are visible here.
task_dir=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(task_dir),'common'));
if ~exist('review_options','var'); review_options=struct(); end
if ~isfield(review_options,'Mode')
    review_options.Mode='reference';
    if isfield(review_options,'Quick') && review_options.Quick; review_options.Mode='compute'; end
end
ctx=review_setup(mfilename('fullpath'),'plot_figure2',review_options);

%% 1. Choose the numerical result explicitly.
if strcmp(ctx.mode,'reference') || strcmp(ctx.mode,'check')
    source_file=fullfile(ctx.data_dir,'reference','result.mat');
else
    source_file=review_find_result(task_dir,'search',ctx.options);
    if isempty(source_file) && ~isfield(review_options,'RunId') && isempty(ctx.options.InputDir)
        source_file=review_find_result(task_dir,'search',rmfield(ctx.options,'RunId'));
    end
    assert(~isempty(source_file),'review:MissingSearch','Run search.m first, or set InputDir to its run.');
end
saved=load(source_file,'route'); route=saved.route;
if strcmp(ctx.mode,'check')
    onsager_check_paper_result(route);
    review_result=review_finish(ctx,struct('status','PASS')); return;
end
if ~route.completed
    save(fullfile(ctx.run_dir,'result.mat'),'source_file');
    review_result=review_finish(ctx,struct('status','SKIPPED', ...
        'reason','The actual search is incomplete; no four-state figure can be drawn.', ...
        'source_file',source_file));
    return;
end
paper=onsager_check_paper_result(route);
grid=load(fullfile(ctx.data_dir,'input','grid.mat'),'xyz');
xyz=grid.xyz;

%% 2. Draw each spherical density and its own density color scale.
visibility='off'; if ctx.visible; visibility='on'; end
fig=figure('Visible',visibility,'Color','white','Units','pixels', ...
    'Position',[100 100 1680 470],'InvertHardcopy','off');
[xx,yy,zz]=sphere(50);
for k=1:4
    node=route.nodes(k); left=0.012+(k-1)*0.25;
    ax=axes('Parent',fig,'Units','normalized','Position',[left .205 .215 .72]);
    scatter3(ax,xyz(:,1),xyz(:,2),xyz(:,3),560,node.logrho,'.'); hold(ax,'on');
    surf(ax,.98*xx,.98*yy,.98*zz,'CData',zeros(size(zz))+min(node.logrho)-.1,'EdgeColor','none');
    shading(ax,'flat'); axis(ax,'equal'); axis(ax,'off'); view(ax,3);
    xlim(ax,[-1.12 1.12]); ylim(ax,[-1.12 1.12]); zlim(ax,[-1.12 1.12]);
    colormap(ax,parula(256));
    cb=colorbar(ax,'Location','eastoutside'); cb.Units='normalized';
    cb.Position=[left+.187 .36 .008 .40]; cb.Color=[0 0 0];
    cb.Ticks=[min(node.logrho),max(node.logrho)];
    cb.TickLabels=compose('%.3g',cb.Ticks); cb.FontSize=11; cb.TickLength=0;
    ax.Position=[left .23 .185 .68];
    annotation(fig,'textbox',[left+.17 .275 .06 .06],'String','log(\rho)', ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',11,'Color','black');
    annotation(fig,'textbox',[left .085 .215 .07],'String',sprintf('index %d',node.index), ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',17,'FontWeight','bold','Color','black');
    annotation(fig,'textbox',[left .01 .215 .08],'String',sprintf('E = %.6f',node.energy), ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',17,'Color','black');
    if k<4
        x=left+.224;
        annotation(fig,'arrow',[x,x+.023],[.57,.57],'LineWidth',1.7,'HeadWidth',8,'HeadLength',8,'Color','black');
    end
end

%% 3. Export the figure and record exactly which route supplied it.
drawnow;
exportgraphics(fig,fullfile(ctx.run_dir,'figure2.png'),'Resolution',220,'BackgroundColor','white','Padding',16);
exportgraphics(fig,fullfile(ctx.run_dir,'figure2.pdf'),'ContentType','image', ...
    'Resolution',300,'BackgroundColor','white','Padding',16);
if ~ctx.visible; close(fig); end

%% 4. Export each state separately with the same spherical density plot.
state_dir=fullfile(ctx.run_dir,'states');
if ~isfolder(state_dir); mkdir(state_dir); end
state_images=cell(1,numel(route.nodes));
for k=1:numel(route.nodes)
    node=route.nodes(k);
    state_fig=figure('Visible',visibility,'Color','white','Units','pixels', ...
        'Position',[100 100 700 650],'InvertHardcopy','off');
    ax=axes('Parent',state_fig,'Units','normalized','Position',[.04 .16 .76 .78]);
    scatter3(ax,xyz(:,1),xyz(:,2),xyz(:,3),560,node.logrho,'.'); hold(ax,'on');
    surf(ax,.98*xx,.98*yy,.98*zz,'CData',zeros(size(zz))+min(node.logrho)-.1,'EdgeColor','none');
    shading(ax,'flat'); axis(ax,'equal'); axis(ax,'off'); view(ax,3);
    xlim(ax,[-1.12 1.12]); ylim(ax,[-1.12 1.12]); zlim(ax,[-1.12 1.12]);
    colormap(ax,parula(256));
    cb=colorbar(ax,'Location','eastoutside'); cb.Color=[0 0 0];
    cb.Ticks=[min(node.logrho),max(node.logrho)];
    cb.TickLabels=compose('%.3g',cb.Ticks); cb.FontSize=13; cb.TickLength=0;
    cb.Label.String='log(\rho)';
    annotation(state_fig,'textbox',[.03 .025 .94 .08], ...
        'String',sprintf('index %d | E = %.6f',node.index,node.energy), ...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',17,'Color','black');
    state_images{k}=fullfile(state_dir,sprintf('index%d.png',node.index));
    drawnow;
    exportgraphics(state_fig,state_images{k},'Resolution',220,'BackgroundColor','white','Padding',16);
    if ~ctx.visible; close(state_fig); end
end
save(fullfile(ctx.run_dir,'result.mat'),'source_file','paper','state_images');
review_result=review_finish(ctx,struct('status','PASS','source_file',source_file, ...
    'indices',[route.nodes.index],'state_images',{state_images},'numerical_search_performed',false));
