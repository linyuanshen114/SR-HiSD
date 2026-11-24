function v=lb_normalize_rms(v,innp)
%LB_NORMALIZE_RMS Original benchmark's zero-mean RMS direction normalization.
% The archived raw HEX state-kick vector deliberately bypasses this helper.
v=v(:); v=v-mean(v); v=v/max(sqrt(innp(v,v)),1e-14);
end
