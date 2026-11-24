function [slice,u,v]=lb_oblique_slice(phi,L,resolution)
%LB_OBLIQUE_SLICE Periodic trilinear physical slice for ordinary state panels.
if nargin<3, resolution=256; end
N=round(numel(phi)^(1/3)); assert(N^3==numel(phi));
e1=[0 1 1]/sqrt(2); e2=[1 1 -1]/sqrt(3);
u=linspace(-L/2,L/2,resolution); v=u; [U,V]=meshgrid(u,v);
X=mod(L/2+e1(1)*U+e2(1)*V,L);
Y=mod(L/2+e1(2)*U+e2(2)*V,L);
Z=mod(L/2+e1(3)*U+e2(3)*V,L);
a=reshape(phi,N,N,N); a=a([1:N 1],[1:N 1],[1:N 1]); grid=(0:N)*L/N;
F=griddedInterpolant({grid,grid,grid},a,'linear','none'); slice=F(X,Y,Z);
end
