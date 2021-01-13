% load('bunny_h0.01.mat');
% load('bunny_full_res.mat');
% load('../fields/tooth8_h0.10.mat')
% load('frog');
% load('happy_buddha.mat');
% load('../fields/stanford_dragon_highres.mat');
% load('../fields/stanford_dragon_fullres_nrml30.mat');
% load('../fields/happy_buddha_highres.mat');
load('../fields/homer.mat');
% load('../fields/armadillo.mat');
%load('../fields/mammoth_tooth.mat');
% load('../fields/cantius_tooth.mat');
% load('../fields/interlocked_tori_big.mat');
% load('../fields/Laurent_Hand.mat');
% load('../fields/gargoyle.mat');
% load('../fields/dancing_children.mat');
% load('../fields/trefoil_N2400.mat');

x = xc;
N = size(x,1);
% % Normal vectors to interpolate for the sphere
% u =@(x) [x(:,1) x(:,2) x(:,3)];
% Normal vectors to interpolate for the tooth
% F = @(x,y,z) x.^8 + y.^8 + z.^8 - (x.^2 + y.^2 + z.^2);
% u = @(x) [2*x(:,1).*(4*x(:,1).^6-1), 2*x(:,2).*(4*x(:,2).^6-1), 2*x(:,3).*(4*x(:,3).^6-1)];
% nrml = u(x);
% nrml = nrml./(sqrt(sum(nrml.^2,2)));
% id = 1:N;
% rng(1378271);
% nrml(id,:) = nrml(id,:) + 1e-1*randn(size(nrml(id,:)));

nrml = nrml./(sqrt(sum(nrml.^2,2)));

% Set up meshgrid for contour plot
[minxx,maxxx] = bounds(x);

x = x - minxx;
x = x./max(maxxx-minxx);

[minx,maxx] = bounds(x);

% mm = 56;
mm = 160;
dx = (maxx-minx)/mm;
[X, Y, Z] = meshgrid(linspace(minx(1)-dx(1),maxx(1)+dx(1),mm),linspace(minx(2)-dx(2),maxx(2)+dx(2),mm),linspace(minx(3)-dx(3),maxx(3)+dx(3),mm));

% Add a point on the boundary of the shape
% for setting the potential constant
xe = [X(:) Y(:) Z(:)];
% xe = load('tooth8_h0.02.mat','xc'); xe = xe.xc;
m = length(xe);

%%
% Smoothing parameter:
lambda = 1e1;
regularization = 0;
exactinterp = 1;

%Define eta and zeta
phi = @(r) r;
% eta = @(r) r.^3;
% zeta = @(r) 3*r;  order = 2;
eta = @(r) -r;
zeta = @(r) -1./r; order = 1;

n = length(x);    % Number of nodes on patch

% Calculate the distance matrix between the nodes on each patch
e = ones(1,n);
dx = (x(:,1) - x(:,1).');
dy = (x(:,2) - x(:,2).');
dz = (x(:,3) - x(:,3).');
r = sqrt(dx.^2 + dy.^2 + dz.^2);

%Vector field coefficient vector
ui = nrml;
b = zeros(3*n,1);
b(1:3:3*n,1) = ui(:,1);
b(2:3:3*n,1) = ui(:,2);
b(3:3:3*n,1) = ui(:,3);

%Construct A
A = zeros(3*n,3*n);
eta_temp = eta(r);
zeta_temp = zeta(r);
zeta_temp(r == 0) = 0;
dphi_xx = zeta_temp.*dx.^2+eta_temp;
dphi_yy = zeta_temp.*dy.^2+eta_temp;
dphi_zz = zeta_temp.*dz.^2+eta_temp;

dphi_xy = zeta_temp.*dx.*dy;
dphi_xz = zeta_temp.*dx.*dz;
dphi_yz = zeta_temp.*dy.*dz;

A(1:3:3*n,1:3:3*n) = dphi_xx;
A(1:3:3*n,2:3:3*n) = dphi_xy;
A(1:3:3*n,3:3:3*n) = dphi_xz;

A(2:3:3*n,1:3:3*n) = dphi_xy;
A(2:3:3*n,2:3:3*n) = dphi_yy;
A(2:3:3*n,3:3:3*n) = dphi_yz;

A(3:3:3*n,1:3:3*n) = dphi_xz;
A(3:3:3*n,2:3:3*n) = dphi_yz;
A(3:3:3*n,3:3:3*n) = dphi_zz;

[CFP,P] = curlfreePoly(x,order);
if regularization == 1    % gcv regularization
    % [U,~,~] = svd(CFP);
    % ZZ = U(:,10:end);
    % C = U(:,1:10)*U(:,1:10)';
    % AZ = A*ZZ;
    % [Q,D] = eig(ZZ'*AZ);
    % D = diag(D);
    % ZQ = ZZ*Q;
    % AZQ = AZ*Q;
    % [lam2,fval,flag,output] = fminbnd(@gcv_cost_phs2,0,35,[],ZQ,AZQ,C,D,b,3*n);
    
    [F1,G1] = qr(CFP);
    F2 = F1(:,10:end);
    F1 = F1(:,1:9);
    G1 = G1(1:9,1:9);
    w1 = F1'*b;
    w2 = F2'*b;
    L = chol(F2'*A*F2);
    [U,D,V] = svd(L');
    D = diag(D);
    z = U'*w2;
    
    % Determine the parameter
    [lam,fval,flag,output] = fminbnd(@gcv_cost_phs3,-10,35,[],z,D,3*n);
    lams(k) = 3*n*exp(-lam);
    A = A + lams(k)*eye(3*n);
    
    % Compute the coefficients
    coeffs = F2*(U*(z./(D.^2 + lams(k))));
    coeffsp = G1\(w1-F1'*(A*coeffs));
else
    if regularization == 2   % Manual regularization
        A = A + lambda*eye(3*n);
    end
    % Compute the coefficients
    l = size(CFP,2);
    coeffs = [A CFP;CFP.' zeros(l)]\[b;zeros(l,1)];
    coeffsp = coeffs(3*n+1:end);
    coeffs = coeffs(1:3*n);
end

temp_potential_nodes = sum(eta_temp.*(dx.*coeffs(1:3:3*n)' + dy.*coeffs(2:3:3*n)' + dz.*coeffs(3:3:3*n)'),2) + P*coeffsp;

if exactinterp
    P = ones(n,1);
    coeffs_correction = [phi(r) P;P.' zeros(1)]\[temp_potential_nodes;zeros(1,1)];
else
    P = [ones(n,1) P(:,1:3)];
    coeffs_correction = P\temp_potential_nodes;
end


% Generate the distance matrix between the patch nodes and global evaluation points
mm = length(xe);

batch_sz = ceil(800^2/n);
temp_potential = zeros(mm,1);
potential_correction = zeros(mm,1);
for j = 1:batch_sz:mm
    idb = j:min((j+batch_sz-1),mm);
    mmb = length(idb);
    dx = (xe(idb,1) - x(:,1)');
    dy = (xe(idb,2) - x(:,2)');
    dz = (xe(idb,3) - x(:,3)');
    r = sqrt(dx.^2 + dy.^2 + dz.^2);
    
    [~,P] = curlfreePoly(xe(idb,:),order);
    
    temp_potential(idb) = sum(eta(r).*(dx.*coeffs(1:3:3*n)' + dy.*coeffs(2:3:3*n)' + dz.*coeffs(3:3:3*n)'),2) + P*coeffsp;
    
    % Potential
    if exactinterp
        potential_correction(idb) = [phi(r) ones(mmb,1)]*coeffs_correction;
    else
        potential_correction(idb) = [ones(mmb,1) P(:,1:3)]*coeffs_correction;
    end
end

% Correct the potential function
g = temp_potential - potential_correction;
g = reshape(g,size(X));

%%
fv = isosurface(X,Y,Z,g,0);
p = patch(fv);
isonormals(X,Y,Z,g,p);
daspect([1 1 1])
view(3)
% set(p,'FaceColor',[255 241 0]/255,'EdgeColor','none')
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
% set(p,'SpecularColorReflectance', 0, 'SpecularExponent', 100, 'DiffuseStrength', 0.8);
lighting phong
material dull
% view([180 -90])
% view([0 90])
% view([60 15])    % Dragon
% view([-90 4]);    % Armadillo
% view([90 0]);   % Buddha
% view([-128 36]);    % cantius_tooth
% camlight headlight
camlight left
axis off
set(gcf,'Color',[1 1 1]*1)
% view([-40 20])
% view([-30 10]);
axis tight
ax = axis;

% %%
% fv = isosurface(X,Y,Z,g,0);
% p = patch(fv);
% isonormals(X,Y,Z,g,p);
% daspect([1 1 1])
% view(3)
% % set(p,'FaceColor',[255 241 0]/255,'EdgeColor','none')
% set(p,'FaceColor',0.95*[1 1 1],'EdgeColor','none')
% % set(p,'SpecularColorReflectance', 0, 'SpecularExponent', 100, 'DiffuseStrength', 0.8);
% material dull
% lighting phong
% % view([180 -90])
% camlight headlight
% axis off
% set(gcf,'Color',[1 1 1]*1)
%
% hold on
% uu = nrml(:,1);
% vv = nrml(:,2);
% ww = nrml(:,3);
% G = sqrt(uu.^2 + vv.^2 + ww.^2);
% % [Fv,Vv,Cv]=quiver3Dpatch(x(:,1),x(:,2),x(:,3),uu,vv,ww, ones(size(x)),0.33*[0 1]);
% [Fv,Vv,Cv]=quiver3Dpatch(x(:,1),x(:,2),x(:,3),uu,vv,ww, ones(size(x)),[0.2 0.2 0.5]);
% patch('Faces',Fv,'Vertices',Vv,'EdgeColor','none', 'CData',Cv,'FaceColor',[1 0 0],'FaceAlpha',1);
% % quiver3(x(:,1),x(:,2),x(:,3),uu,vv,ww,1,'k-')
% % plot3(x(:,1),x(:,2),x(:,3),'k.')
% daspect([1 1 1])
% % axis(ax)
% view(3)
% axis off

