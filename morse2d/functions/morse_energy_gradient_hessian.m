function [V, g, H] = morse_energy_gradient_hessian(q, a, r0)
%MORSE_ENERGY_GRADIENT_HESSIAN Pair energy, gradient and exact Cartesian Hessian.
    if nargin < 3, r0 = 1; end
    if nargin < 2, a = 6;  end

    nd = length(q);
    if mod(nd,2) ~= 0, error('q length must be multiple of 2'); end
    N = nd/2;
    pos = reshape(q,2,N)';

    V = 0;
    g = zeros(nd,1);
    H = zeros(nd,nd);
    I2 = eye(2);

    for i=1:N-1
        ii=(i-1)*2+1:i*2; pi=pos(i,:)';
        for j=i+1:N
            jj=(j-1)*2+1:j*2; pj=pos(j,:)';
            d = pi-pj;
            r2 = d.'*d;
            r = sqrt(r2);
            if r < 1e-12, r=1e-12; r2=r*r; end

            t  = exp(-a*(r-r0));
            t2 = t*t;
            V = V + (t2 - 2*t);

            dVdr = 2*a*(t - t2);
            gs = dVdr/r;
            gij = gs*d;
            g(ii)=g(ii)+gij;
            g(jj)=g(jj)-gij;

            if nargout>2
                d2 = 2*a*a*(2*t2 - t);
                A = (d2 - gs)/r2;
                B = gs;
                Hp = A*(d*d.') + B*I2;
                H(ii,ii)=H(ii,ii)+Hp;
                H(jj,jj)=H(jj,jj)+Hp;
                H(ii,jj)=H(ii,jj)-Hp;
                H(jj,ii)=H(jj,ii)-Hp;
            end
        end
    end
end
