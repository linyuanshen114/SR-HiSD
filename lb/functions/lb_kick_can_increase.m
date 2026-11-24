function b=lb_kick_can_increase(force_below_tol,gate_verified,candidate_index,parent_index,field_rms,kick,energy,parent_energy)
%LB_KICK_CAN_INCREASE Permit a larger kick only for an ineffective branch.
b=force_below_tol || (gate_verified && candidate_index==parent_index && field_rms<=.25*kick ...
    && abs(energy-parent_energy)<=max(1e-10,1e-9*abs(parent_energy)));
end
