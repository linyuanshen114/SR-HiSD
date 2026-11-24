function info=lb_branch_info(phi,bcc,parent,sign,hex_reference)
%LB_BRANCH_INFO Residual and conservative phase label for one actual GF end.
if nargin<5; hex_reference=[]; end
force=lb_negative_gradient(phi); rms=sqrt(mean(force.^2));
N=round(numel(phi)^(1/3)); a=reshape(phi,N,N,N); b=reshape(bcc,N,N,N);
cross=real(ifftn(fftn(b).*conj(fftn(a)))); [~,id]=max(cross(:));
[i,j,k]=ind2sub([N N N],id); aligned=circshift(a,[i j k]-1);
distance=sqrt(mean((aligned(:)-bcc(:)).^2));
label='unconverged';
if rms<1e-9
    label='Other endpoint';
    if distance<1e-3 && abs(lb_energy(phi)-lb_energy(bcc))<1e-8, label='BCC';
    elseif sqrt(mean(phi.^2))<1e-4 && abs(lb_energy(phi))<1e-9, label='disorder'; end
end
phase_evidence=struct('hex_phase_certified',false);
if strcmp(label,'Other endpoint') && ~isempty(hex_reference)
    phase_evidence=lb_hex_identity(phi,hex_reference.phi);
    if phase_evidence.hex_phase_certified, label='HEX'; end
end
info=struct('phi',phi,'sign',sign,'class',label,'residual',rms,'converged',rms<1e-9, ...
    'parent_phi_sha256',lb_numeric_sha(parent),'phi_sha256',lb_numeric_sha(phi), ...
    'phase_evidence',phase_evidence);
end
