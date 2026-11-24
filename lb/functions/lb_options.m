function options=lb_options(k,V,bb,N)
%LB_OPTIONS Translate the explicit BB/settings structure into kernel options.
options=lb_solver_options(1);
options.k=k; options.V=V; options.dt=bb.dt;
options.betau=bb.betau; options.betat=bb.betat;
options.maxiter=bb.maxiter; options.epsf=bb.epsf;
options.outputp=bb.outputp; options.outputd=Inf; options.drawfun=[];
options.innpfunc=@(X,Y) (X'*Y)/N^3; options.precd=@lb_precondition;
end
