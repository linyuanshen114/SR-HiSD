%% Combine the seven Table 1 rows and draw the success-rate matrix
% In compute mode, input rates must be computed in the same run. Reference
% mode loads the current saved full-calculation rates.
if ~exist('review_options','var'), review_options=struct(); end
task_dir=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(task_dir),'common'));
plot_options=review_options;
if ~isfield(plot_options,'Mode')
    plot_options.Mode='reference';
    if isfield(plot_options,'Quick') && plot_options.Quick; plot_options.Mode='compute'; end
end
ctx=review_setup(mfilename('fullpath'),'plot_table1',plot_options);
methods=["1-HiSD";"2-HiSD";"3-HiSD";"4-HiSD"; ...
    "NPSS, epsilon_lambda=0.10";"NPSS, epsilon_lambda=0.01";"1-SR-HiSD"];
run_info=struct('subtask','plot_table1','mode',ctx.mode,'quick',ctx.quick, ...
    'actual_solver_runs',0,'methods',methods);
if strcmp(ctx.mode,'check')
    assert(isfile(fullfile(ctx.data_dir,'paper','table1.csv')));
    review_result=review_finish(ctx,run_info); return;
end

%% Load rates from this calculation or the current saved full calculation
if strcmp(ctx.mode,'compute')
    allow_latest=~isfield(review_options,'RunId');
    [original_file,npss_file]=morse_table_sources(ctx,allow_latest);
    original=load(original_file,'Rates','run_info'); npss=load(npss_file,'Rates','run_info');
    Rates=[original.Rates;npss.Rates(:,original.Rates.Properties.VariableNames)];
    Ns=original.run_info.Ns;
    source_files={original_file,npss_file};
else
    original_file=fullfile(ctx.data_dir,'reference','original_methods','result.mat');
    npss_file=fullfile(ctx.data_dir,'reference','npss','result.mat');
    original=load(original_file,'Rates','run_info'); npss=load(npss_file,'Rates','run_info');
    Rates=[original.Rates;npss.Rates(:,original.Rates.Properties.VariableNames)];
    Ns=original.run_info.Ns;
    source_files={original_file,npss_file};
end

%% Build the table in the paper's explicit method order
percentages=nan(numel(methods),numel(Ns));
for method_id=1:numel(methods)
    for column=1:numel(Ns)
        row=Rates(Rates.method==methods(method_id) & Rates.N==Ns(column),:);
        assert(height(row)==1 && row.trials>0);
        percentages(method_id,column)=100*row.successes/row.trials;
    end
end
Table1=array2table(percentages,'VariableNames',cellstr("N"+string(Ns)));
Table1=addvars(Table1,methods,'Before',1,'NewVariableNames','method');
if ~ctx.quick
    expected=readtable(fullfile(ctx.data_dir,'paper','table1.csv'),'TextType','string');
    assert(isequal(Table1.method,expected.method) && ...
        max(abs(Table1{:,2:end}-expected{:,2:end}),[],'all')<1e-10, ...
        'Review:MorseTableMismatch','Computed rates differ from the paper Table 1; actual measurements are retained.');
end

%% Plot the measured percentages; labels and values are drawn here
fig=figure('Color','w','Visible',on_off(ctx.visible),'Position',[100 100 1050 520]);
ax=axes(fig,'Position',[.31 .15 .60 .73]);
imagesc(ax,1:numel(Ns),1:numel(methods),percentages,[0 100]);
colormap(ax,parula(256));
ax.YDir='reverse'; ax.XTick=1:numel(Ns); ax.XTickLabel="N="+string(Ns);
ax.YTick=1:numel(methods); ax.YTickLabel=methods; ax.TickLabelInterpreter='none';
ax.FontSize=12; ax.XColor=[.1 .1 .1]; ax.YColor=[.1 .1 .1];
for row=1:numel(methods)
    for column=1:numel(Ns)
        value=percentages(row,column); color=[1 1 1];
        if value>60, color=[.1 .1 .1]; end
        text(ax,column,row,sprintf('%.1f',value),'HorizontalAlignment','center', ...
            'FontSize',12,'FontWeight','bold','Color',color);
    end
end
bar=colorbar(ax); bar.Label.String='Index-one success (%)'; bar.Color=[.1 .1 .1];
label='Table 1: Morse clusters, shared initial conditions';
if ctx.quick, label=[label ' (bounded execution check)']; end
title(ax,label,'Color',[.1 .1 .1],'FontSize',14);
exportgraphics(fig,fullfile(ctx.run_dir,'table1.png'),'Resolution',220);
exportgraphics(fig,fullfile(ctx.run_dir,'table1.pdf'),'ContentType','vector');
if ~ctx.visible, close(fig); end

%% Save numerical table and its source association separately from the image
writetable(Table1,fullfile(ctx.run_dir,'table1.csv'));
run_info.Ns=Ns; run_info.source_files=source_files;
save(fullfile(ctx.run_dir,'result.mat'),'Table1','Rates','percentages','run_info','-v7');
disp(Table1);
review_result=review_finish(ctx,run_info);

function value=on_off(visible)
if visible, value='on'; else, value='off'; end
end
