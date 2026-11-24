function selected=lb_choose_candidate(candidates,target)
%LB_CHOOSE_CANDIDATE First accepted sign; otherwise retain a labeled failure.
for j=1:numel(candidates)
    if candidates{j}.valid, selected=candidates{j}; return; end
end
score=cellfun(@(c) 1e6*abs(c.index-target)+c.residual,candidates);
[~,j]=min(score); selected=candidates{j};
end
