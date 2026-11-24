function [f,raw,energy] = onsager_force(logrho,U,p)
rho = exp(logrho);
Krho = p.Kw*rho;
raw = -logrho-Krho;
raw = raw-mean(raw);
f = onsager_project(raw,U,p.w);
energy = p.w*(rho'*logrho+0.5*rho'*Krho);
end

