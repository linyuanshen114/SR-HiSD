function evidence=lb_hex_identity(phi,reference_phi)
%LB_HEX_IDENTITY Numerical HEX phase identity against an independent field.
% The N=128 morphology uses the six-wave triangular star, phase-locked
% triad, cylindrical invariance and represented-mode sixfold symmetry.
% Proper cubic FFT index rotations and a continuous periodic translation
% are allowed; field relaxation and amplitude rescaling are not performed.
% This phase test is not a Hessian/inertia or minimum-stability certificate.
N=128;
assert(isnumeric(phi) && isreal(phi) && numel(phi)==N^3 && all(isfinite(phi(:))), ...
    'LB:HEXField','HEX identity requires a finite real 128^3 field.');
assert(isnumeric(reference_phi) && isreal(reference_phi) && numel(reference_phi)==N^3 && ...
    all(isfinite(reference_phi(:))),'LB:HEXReference','HEX reference must be finite, real and 128^3.');
phi=double(phi(:)); reference_phi=double(reference_phi(:));
proper_reference_rotation=[0 1 0;1 0 0;0 0 -1];
empty_alignment=struct('continuous_grid_shift',[NaN NaN NaN],'rms_difference',Inf, ...
    'relative_rms_difference',Inf,'method', ...
    'two independent primary-wave phases; full-spectrum Parseval residual, no amplitude rescale', ...
    'global_optimum_claim',false,'source_or_target_modified',false);
evidence=struct('phase','unresolved','numerical_phase','unresolved', ...
    'hex_phase_certified',false,'spectrum',empty_geometry(), ...
    'reference_alignment',empty_alignment,'endpoint_proper_rotation',eye(3), ...
    'proper_reference_rotation',proper_reference_rotation, ...
    'reference_relative_rms_tolerance',1e-5,'reference_morphology_pass',false, ...
    'phi_sha256',lb_numeric_sha(phi),'reference_phi_sha256',lb_numeric_sha(reference_phi), ...
    'strict_minimum_hessian_certified',false,'source_field_modified',false, ...
    'field_relaxation_performed',false,'reason','');
p=reshape(reference_phi,N,N,N);
oriented=permute(p,[2 1 3]); oriented=oriented(:,:,[1 N:-1:2]);
Uref=fftn(oriented-mean(oriented(:)))/(N^3);
power=abs(Uref).^2; power(1)=0; total=sum(power(:));
if ~isfinite(total) || total<=0
    evidence.reason='Independent reference has no nonconstant spectral power.'; return;
end
reference_check=hex_geometry(Uref);
evidence.reference_morphology_pass=reference_check.morphology_pass;
evidence.reference_spectrum=reference_check;
if ~reference_check.morphology_pass
    evidence.reason='Independent reference failed HEX geometry.'; return;
end
v=reshape(phi,N,N,N); U0=fftn(v-mean(v(:)))/(N^3);
power=abs(U0).^2; power(1)=0; total=sum(power(:));
if ~isfinite(total) || total<=0
    evidence.reason='Endpoint has no nonconstant spectral power.'; return;
end
[U,shape,endpoint_rotation]=hex_orient(U0);
alignment=hex_align(Uref,U);
matched=alignment.relative_rms_difference<1e-5;
passed=shape.morphology_pass && matched;
evidence.spectrum=shape; evidence.reference_alignment=alignment;
evidence.endpoint_proper_rotation=endpoint_rotation;
evidence.hex_phase_certified=logical(passed);
if passed
    evidence.phase='HEX'; evidence.numerical_phase='HEX'; evidence.reason='';
elseif ~shape.morphology_pass
    evidence.reason='Endpoint failed HEX geometry.';
else
    evidence.reason='Endpoint did not match the independent HEX reference.';
end
end

function s=empty_geometry()
checks=struct('principal_star_matches',false,'principal_rank2',false, ...
    'equal_radius',false,'primary_power',false,'amplitude_equality',false, ...
    'cylindrical_axis_invariance',false,'triangular_reciprocal_lattice',false, ...
    'positive_cylinder_triad',false,'physical_hex_sixfold',false,'sixfold_coverage',false);
s=struct('primary_vectors',zeros(0,3),'primary_amplitudes',zeros(0,1), ...
    'column_axis',[1 1 -1],'primary_power_fraction',0, ...
    'primary_amplitude_relative_spread',NaN,'off_plane_power_fraction',NaN, ...
    'off_lattice_power_fraction',NaN,'triad_phase_error_radians',NaN, ...
    'sixfold_relative_defect',Inf,'sixfold_covered_power_fraction',0, ...
    'checks',checks,'morphology_pass',false);
end

function s=hex_geometry(U)
N=size(U,1); q1=[12 -12 0]; q2=[0 12 12]; expected=[q1;-q1;q2;-q2;q1+q2;-q1-q2];
f=[0:N/2-1 -N/2:-1]; power=abs(U).^2; power(1)=0; total=sum(power(:));
assert(isfinite(total) && total>0,'No nonconstant spectral power.');
[~,ids]=maxk(power(:),6); [a,b,c]=ind2sub([N N N],ids);
vectors=sortrows([reshape(f(a),[],1) reshape(f(b),[],1) reshape(f(c),[],1)]);
amplitudes=abs(U(hex_kindex(vectors,N))); k2=sum(vectors.^2,2);
same_star=isequal(vectors,sortrows(expected));
[I,J,K]=ndgrid(f,f,f); plane=I+J-K==0;
lattice=plane & mod(I,12)==0 & mod(J,12)==0 & mod(K,12)==0;
off_plane=sum(power(~plane))/total; off_lattice=sum(power(~lattice))/total;
primary=sum(abs(U(hex_kindex(expected,N))).^2)/total;
spread=(max(amplitudes)-min(amplitudes))/mean(amplitudes);
p1=angle(U(hex_kindex(q1,N))); p2=angle(U(hex_kindex(q2,N)));
triad=abs(angle(exp(1i*(p1+p2-angle(U(hex_kindex(q1+q2,N)))))));
[a,b,c]=ind2sub([N N N],find(lattice));
ks=[reshape(f(a),[],1) reshape(f(b),[],1) reshape(f(c),[],1)];
m=ks(:,1)/12; n=ks(:,3)/12; rotated=(m-n)*q1+m*q2;
safe=all(rotated>=-N/2 & rotated<N/2,2);
ks=ks(safe,:); rotated=rotated(safe,:); m=m(safe); n=n(safe);
A=U(hex_kindex(ks,N)).*exp(-1i*(m*p1+n*p2));
B=U(hex_kindex(rotated,N)).*exp(-1i*(rotated(:,1)/12*p1+rotated(:,3)/12*p2));
covered=sum(abs(A).^2)/total; defect=norm(A-B)/norm(A);
checks=struct('principal_star_matches',same_star,'principal_rank2',rank(vectors)==2, ...
    'equal_radius',all(k2==288),'primary_power',primary>.99, ...
    'amplitude_equality',spread<1e-5,'cylindrical_axis_invariance',off_plane<1e-10, ...
    'triangular_reciprocal_lattice',off_lattice<1e-10, ...
    'positive_cylinder_triad',triad<1e-6,'physical_hex_sixfold',defect<1e-5, ...
    'sixfold_coverage',covered>=1-1e-10);
s=struct('primary_vectors',vectors,'primary_amplitudes',amplitudes,'column_axis',[1 1 -1], ...
    'primary_power_fraction',primary,'primary_amplitude_relative_spread',spread, ...
    'off_plane_power_fraction',off_plane,'off_lattice_power_fraction',off_lattice, ...
    'triad_phase_error_radians',triad,'sixfold_relative_defect',defect, ...
    'sixfold_covered_power_fraction',covered,'checks',checks, ...
    'morphology_pass',all(structfun(@(x) logical(x),checks)));
end

function info=hex_align(source,target)
N=size(source,1); A=[12 -12 0;0 12 12];
phase=angle(target(hex_kindex(A,N))./source(hex_kindex(A,N)));
shift=pinv(A)*(phase*N/(2*pi));
f=[0:N/2-1 -N/2:-1];
phase_grid=(2*pi/N)*(reshape(f,[],1,1)*shift(1)+reshape(f,1,[],1)*shift(2)+reshape(f,1,1,[])*shift(3));
difference=source.*exp(1i*phase_grid)-target;
rms=norm(difference(:));
info=struct('continuous_grid_shift',shift(:).','rms_difference',rms, ...
    'relative_rms_difference',rms/norm(target(:)), ...
    'method','two independent primary-wave phases; full-spectrum Parseval residual, no amplitude rescale', ...
    'global_optimum_claim',false,'source_or_target_modified',false);
end

function [U,shape,R]=hex_orient(original)
% First keep the archived orientation untouched. Fresh runs can select an
% equivalent cubic orientation; compare them in one canonical frame using
% only unitary FFT index permutations, never by interpolating/relaxing phi.
R=eye(3); U=original; shape=hex_geometry(U);
if ~shape.morphology_pass
    candidates={[0 1 0;1 0 0;0 0 -1]};
    permutations=perms(1:3);
    for i=1:size(permutations,1)
        p=permutations(i,:);
        for bits=0:7
            signs=2*bitget(bits,1:3)-1;
            Q=zeros(3); Q(sub2ind([3 3],(1:3)',p(:)))=signs(:);
            if round(det(Q))==1; candidates{end+1}=Q; end %#ok<AGROW>
        end
    end
    for i=1:numel(candidates)
        trial=hex_rotate_spectrum(original,candidates{i});
        if ~hex_star_matches(trial); continue; end
        check=hex_geometry(trial);
        if check.morphology_pass
            U=trial; shape=check; R=candidates{i}; break;
        end
    end
end
axis=(R*[1;1;-1]).'; first=find(axis~=0,1);
if axis(first)<0; axis=-axis; end
shape.original_column_axis=axis;
shape.canonical_comparison_axis=[1 1 -1];
end

function yes=hex_star_matches(U)
N=size(U,1); f=[0:N/2-1 -N/2:-1]; p=abs(U).^2; p(1)=0;
[~,ids]=maxk(p(:),6); [a,b,c]=ind2sub([N N N],ids);
vectors=[reshape(f(a),[],1) reshape(f(b),[],1) reshape(f(c),[],1)];
q1=[12 -12 0]; q2=[0 12 12];
yes=isequal(sortrows(vectors),sortrows([q1;-q1;q2;-q2;q1+q2;-q1-q2]));
end

function rotated=hex_rotate_spectrum(U,R)
[~,p]=max(abs(R),[],2); inverse=zeros(1,3); inverse(p)=1:3;
rotated=permute(U,inverse); N=size(U,1);
for axis=1:3
    if R(inverse(axis),axis)<0
        subs={':',':',':'}; subs{axis}=[1 N:-1:2]; rotated=rotated(subs{:});
    end
end
end

function ids=hex_kindex(k,N)
assert(size(k,2)==3); k=mod(k,N)+1;
ids=sub2ind([N N N],k(:,1),k(:,2),k(:,3));
end

