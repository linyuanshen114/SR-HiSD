function k=lb_next_lbsd_k(actual_index)
%LB_NEXT_LBSD_K Drop one direction from the verified actual saddle index.
assert(actual_index>=2 && actual_index<=5 && actual_index==floor(actual_index));
k=actual_index-1;
end
