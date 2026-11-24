function style=lb_plot_style(case_id)
%LB_PLOT_STYLE Final LB phase colors and ordinary linear display ranges.
% The nonuniform color-node positions are part of the published style.
% Values outside CLim saturate at the endpoint colors; data are unchanged.
if nargin<1 || isempty(case_id); case_id='bcc2hex'; end
case_id=lower(char(string(case_id)));
positions=[0 .45 .55 1];
rgb=[.04 .10 .36; .35 .78 1.00; .25 .68 .18; 1.00 .94 .05];
colors=struct('LBSD',[0 117 171]/255,'SRHiSD',[189 53 74]/255, ...
    'GF',[113 119 123]/255,'text',[32 43 50]/255,'guide',[170 179 184]/255);
if contains(case_id,'disorder')
    limits=[-.4 .8]; ticks=[-.4 0 .4 .8]; family='bcc2disorder';
else
    limits=[-1 1]; ticks=[-1 0 1]; family='bcc2hex';
end
style=struct('name','lb_phase','case_id',family,'positions',positions, ...
    'rgb',rgb,'cmap',interp1(positions,rgb,linspace(0,1,256),'linear'), ...
    'clim',limits,'ticks',ticks,'colors',colors, ...
    'normalization','linear','clip',true);
end
