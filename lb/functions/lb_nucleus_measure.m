function measured=lb_nucleus_measure(phi,bcc,other,L,resolution)
%LB_NUCLEUS_MEASURE Port of the final manual endpoint-phase/PCA recipe.
% Inputs are actual same-parent fields. This routine does not classify them.
% It preserves their values and performs display-only mean removal/alignment.
if nargin<5; resolution=1400; end
N=round(numel(phi)^(1/3));
assert(N^3==numel(phi) && numel(bcc)==N^3 && numel(other)==N^3);
assert(all(isfinite([phi(:);bcc(:);other(:)])) && isreal(phi) && isreal(bcc) && isreal(other));
assert(isscalar(L) && L>0 && isfinite(L));
assert(isscalar(resolution) && resolution>=64 && resolution==round(resolution));
means=[mean(phi(:)),mean(bcc(:)),mean(other(:))];
phi=reshape(double(phi(:))-means(1),N,N,N);
bcc=reshape(double(bcc(:))-means(2),N,N,N);
other=reshape(double(other(:))-means(3),N,N,N);
[bcc,bcc_shift]=align_periodic(bcc,phi);
[other,other_shift]=align_periodic(other,phi);
center=max_index_c_order(abs(phi-bcc));
display_shift=[floor(N/2)+1 floor(N/2)+1 floor(N/2)+1]-center;
R=2*sqrt(2)*pi;
axis_distance=min(0:N-1,N-(0:N-1))*L/N;
[X,Y,Z]=ndgrid(axis_distance,axis_distance,axis_distance);
r2=X.^2+Y.^2+Z.^2; kernel=zeros(N,N,N); inside=r2<R^2;
kernel(inside)=exp(R^2./(r2(inside)-R^2));
kernel=kernel/sum(kernel(:)); spectrum=fftn(kernel);
dm=real(ifftn(fftn((phi-bcc).^2).*spectrum));
ds=real(ifftn(fftn((phi-other).^2).*spectrum));
rho=dm./max(dm+ds,eps); rho(~isfinite(rho))=NaN;
phi=circshift(phi,display_shift); rho=circshift(rho,display_shift);
[slice_phi,u,v]=manual_slice(phi,L,resolution);
[slice_rho,~,~]=manual_slice(rho,L,resolution);
[U,V]=meshgrid(u,v);
focus_radius=min(max(u)-min(u),max(v)-min(v))/6;
threshold=.001;
mask=(U.^2+V.^2<=focus_radius^2) & isfinite(slice_rho) & slice_rho>threshold;
ellipse=fit_ellipse(slice_rho,U,V,mask,threshold);
measured=struct('N',N,'L',L,'phi_slice',slice_phi,'phase_density_slice',slice_rho, ...
    'mask',mask,'u',u,'v',v,'ellipse',ellipse,'threshold',threshold, ...
    'mollifier_R',R,'focus_radius',focus_radius,'crop_fraction',.92, ...
    'plane_v1',[0 1 1],'plane_v2',[1 1 -1], ...
    'bcc_alignment_shift',bcc_shift,'other_alignment_shift',other_shift, ...
    'center_index_matlab',center,'center_index_zero_based',center-1, ...
    'display_shift',display_shift,'display_means_removed',means, ...
    'rho_range',[min(rho(:)) max(rho(:))], ...
    'interpolation','periodic cardinal cubic B-spline, matching scipy order=3 grid-wrap', ...
    'method','periodic 3D endpoint phase density and weighted PCA ellipse', ...
    'radius_certificate',false,'source_fields_modified',false);
end

function [aligned,shift]=align_periodic(field,target)
N=size(field,1);
if max(field(:))==min(field(:)) || max(target(:))==min(target(:))
    shift=[0 0 0]; aligned=field; return;
end
correlation=real(ifftn(fftn(target).*conj(fftn(field))));
shift=max_index_c_order(correlation)-1;
shift(shift>N/2)=shift(shift>N/2)-N;
aligned=circshift(field,shift);
end

function sub=max_index_c_order(a)
% NumPy's argmax resolves exact ties in C order; preserve that convention.
reverse=permute(a,[3 2 1]); [~,index]=max(reverse(:));
[k,j,i]=ind2sub(size(reverse),index); sub=[i j k];
end

function [image,u,v]=manual_slice(field,L,resolution)
N=size(field,1); v1=[0 1 1]; v2=[1 1 -1];
e1=v1/norm(v1); v2p=v2-dot(e1,v2)*e1; e2=v2p/norm(v2p);
half_u=.5*L*norm(v1); half_v=.5*L*norm(v2p);
u=-half_u+(0:resolution-1)*(2*half_u/resolution);
v=-half_v+(0:resolution-1)*(2*half_v/resolution);
[U,V]=meshgrid(u,v);
x=mod((L/2+e1(1)*U+e2(1)*V)*N/L,N);
y=mod((L/2+e1(2)*U+e2(2)*V)*N/L,N);
z=mod((L/2+e1(3)*U+e2(3)*V)*N/L,N);
image=periodic_cubic(field,x,y,z);
trim=round((1-.92)*resolution/2);
keep=trim+1:resolution-trim;
image=image(keep,keep); u=u(keep); v=v(keep);
end

function values=periodic_cubic(field,x,y,z)
% Periodic spline prefilter: samples = coefficients * [1,4,1]/6 per axis.
N=size(field,1); h=(2+cos(2*pi*(0:N-1)/N))/3;
denominator=reshape(h,[],1,1).*reshape(h,1,[],1).*reshape(h,1,1,[]);
coeff=real(ifftn(fftn(field)./denominator));
ix=floor(x); iy=floor(y); iz=floor(z);
wx=weights(x-ix); wy=weights(y-iy); wz=weights(z-iz);
values=zeros(size(x));
for a=1:4
    ai=mod(ix+a-2,N)+1;
    for b=1:4
        bj=mod(iy+b-2,N)+1;
        for c=1:4
            ck=mod(iz+c-2,N)+1;
            indices=ai+N*(bj-1)+N^2*(ck-1);
            values=values+coeff(indices).*wx{a}.*wy{b}.*wz{c};
        end
    end
end
end

function w=weights(t)
w={(1-t).^3/6,(3*t.^3-6*t.^2+4)/6, ...
    (-3*t.^3+3*t.^2+3*t+1)/6,t.^3/6};
end

function ellipse=fit_ellipse(rho,U,V,mask,threshold)
ellipse=struct([]);
if nnz(mask)<4; return; end
w=rho(mask)-threshold; x=U(mask); y=V(mask);
good=isfinite(w) & w>0; w=w(good); x=x(good); y=y(good);
if numel(w)<4 || sum(w)<=0; return; end
center=[sum(w.*x) sum(w.*y)]/sum(w); dx=x-center(1); dy=y-center(2);
covariance=[sum(w.*dx.^2) sum(w.*dx.*dy);sum(w.*dx.*dy) sum(w.*dy.^2)]/sum(w);
[vectors,values]=eig((covariance+covariance')/2);
[values,order]=sort(max(real(diag(values)),0),'descend'); vectors=vectors(:,order);
if values(1)<=0; return; end
ellipse=struct('cx',center(1),'cy',center(2),'width',4*sqrt(values(1)), ...
    'height',4*sqrt(values(2)),'angle',atan2(vectors(2,1),vectors(1,1)), ...
    'threshold',threshold,'max_phase_density',max(rho(:)),'mask_points',nnz(mask));
end
