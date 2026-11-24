function HV = onsager_hessian_action(logrho,V,U,p)
HV = V./max(exp(logrho),1e-12)+p.Kw*V;
HV = onsager_project(HV,U,p.w);
end

