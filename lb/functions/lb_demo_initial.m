function seed=lb_demo_initial(case_id,method,kick)
%LB_DEMO_INITIAL Deterministic tiny-grid field for Quick execution checks.
% A bounded test, not a downsampled publication state or convergence claim.
N=8; L=2*pi;
par=struct('tau',-.3,'gamma',.8,'q0',1,'xi',1);
if strcmp(case_id,'bcc2disorder'), par.tau=.05; par.gamma=.6; end
lb_initialize(N,L,par); global LB
bcc=.15*(cos(LB.X).*cos(LB.Y)+cos(LB.X).*cos(LB.Z)+cos(LB.Y).*cos(LB.Z));
bcc=bcc(:)-mean(bcc(:)); rng(17,'twister'); V=randn(N^3,4); V=V-mean(V,1);
[Q,~]=qr(V,0); V=Q*sqrt(N^3);
projection=lb_translation_projector(bcc);
v=projection.PN(V(:,4)); v=v-mean(v); v=v/sqrt(mean(v.^2));
seed=struct('N',N,'L',L,'par',par,'phi_bcc',bcc,'V_bcc_orig',V,'V_bcc_tr',v, ...
    'phi_init',bcc+kick*v,'kick',kick,'down_kick',.005, ...
    'bb',struct('dt',.002,'betau',1,'betat',1,'maxiter',4,'epsf',1e-14,'outputp',4));
if strcmp(method,'lbsd4_to1'), seed.V_bcc_orig(:,4)=v; end
end
