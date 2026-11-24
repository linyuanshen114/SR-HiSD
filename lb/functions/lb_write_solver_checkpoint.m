function lb_write_solver_checkpoint(output_dir,method,state,history)
%LB_WRITE_SOLVER_CHECKPOINT Save this run's current SR state and diagnostics.
global LB
N=LB.N; L=LB.L; par=LB.par; %#ok<NASGU>
temporary=fullfile(output_dir,[method '_checkpoint.tmp.mat']);
destination=fullfile(output_dir,[method '_checkpoint.mat']);
save(temporary,'state','history','N','L','par','-v7.3');
movefile(temporary,destination,'f');
end
