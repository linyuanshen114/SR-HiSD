function [q,res,iterations] = morse_gradient_descent(q,dt,maxiter,a,r0,tol)
%MORSE_GRADIENT_DESCENT Relax the stored minimum with the original Euler step.
% The stopping residual belongs to the gradient before the last Euler step;
% this preserves the archived minimum and the exact shared perturbations.
res=1; iterations=0;
while res>tol && iterations<maxiter
    iterations=iterations+1;
    [~,g]=morse_energy_gradient_hessian(q,a,r0);
    q=q-dt*g;
    res=norm(g)/sqrt(numel(q));
end
end
