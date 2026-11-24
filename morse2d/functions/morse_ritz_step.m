function [vnew,vprevious] = morse_ritz_step(A,v,vprevious)
%MORSE_RITZ_STEP Three-block [direction,residual,previous direction] Ritz update.
v = v/norm(v);
u = A*v;
w = u-v*(v'*u);
w = w-v*(v'*w);
B = v;
AB = u;
if norm(w)>1e-8
    w = w/norm(w);
    B = [B,w];
    AB = [AB,A*w];
end
if ~isempty(vprevious)
    p = vprevious-B*(B'*vprevious);
    if norm(p)>1e-8
        p = p/norm(p);
        p = p-B*(B'*p);
        p = p/norm(p);
        B = [B,p];
        AB = [AB,A*p];
    end
end
K = B'*AB;
[Z,ev] = eig((K+K')/2,'vector');
[~,i] = min(ev);
vnew = B*Z(:,i);
vprevious = v;
end
