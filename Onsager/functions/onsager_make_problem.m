function p = onsager_make_problem(cfg)
S = load(fullfile(cfg.data_dir,'grid.mat'),'xyz');
p.xyz = S.xyz; p.N = size(S.xyz,1); p.w = 4*pi/p.N;
dots = p.xyz*p.xyz';
K = 80*sqrt(abs(1-dots.^2));
% This is the original centering, including its constant energy offset.
K = K-mean(K,2)-mean(K,1)+2*mean(K,'all');
p.Kw = K*p.w;
p.Y = onsager_spherical_matrix(p.xyz,cfg.lmax);
G = p.Y'*(p.Y*p.w);
p.L = chol((G+G')/2+1e-12*eye(size(G)),'lower');
p.Ym = cell(3,1); p.Yp = cell(3,1);
axes = eye(3);
for j = 1:3
    Rm = onsager_rotation_matrix(axes(j,:),-cfg.epsrot);
    Rp = onsager_rotation_matrix(axes(j,:),cfg.epsrot);
    p.Ym{j} = onsager_spherical_matrix((Rm*p.xyz')',cfg.lmax);
    p.Yp{j} = onsager_spherical_matrix((Rp*p.xyz')',cfg.lmax);
end
p.epsrot = cfg.epsrot; p.sr_tol = cfg.sr_tol;
end

