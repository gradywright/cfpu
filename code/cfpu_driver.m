% ptcloud = 'trefoil'
% ptcloud = 'buddha'
% ptcloud = 'homer'
% ptcloud = 'armadillo'
% ptcloud = 'children'
ptcloud = 'dragon'
% ptcloud = 'raptor_head'
% ptcloud = 'filigree'
% ptcloud = 'pump'
% ptcloud = 'gargoyle'
% ptcloud = 'bunny'
% ptcloud = 'hand'

switch ptcloud
    case 'trefoil'
        load('../ptclouds/trefoil_N23064.mat');
        % load('../ptclouds/trefoil_N6144.mat');
        load('../ptclouds/trefoil_patches_M864.mat');
        Q = [0 0 1;1 0 0;0 1 0];
        vw = [90 0];
    case 'buddha'
%         load('../ptclouds/happy_buddha.mat');
        load('../ptclouds/happy_buddha_highres.mat');
%         load('../ptclouds/happy_buddha_highres_patches.mat');
%         load('../ptclouds/happy_buddha_highres_patches_M42861.mat');
        [~,d] = knnsearch(xc,xc,'k',2);
        surface_area = length(xc)*mean(d(:,2))^2;
        y = util.PcCoarsen3D(xc,42861,surface_area);
        Q = [0 0 1;1 0 0;0 1 0];
        vw = [90 15];
    case 'homer'
        load('../ptclouds/homer.mat');
        load('../ptclouds/homer_patches.mat');
        Q = [0 0 1;1 0 0;0 1 0];
        vw = [90 0];
    case 'raptor_head'
        load('../ptclouds/raptor_head.mat');
%         load('../ptclouds/raptor_head_patches.mat');
        load('../ptclouds/raptor_head_patches_M10337.mat');
        Q = [0 0 1;1 0 0;0 1 0];
    case 'raptor'
        load('../ptclouds/raptor.mat');
        load('../ptclouds/raptor_patches.mat');
        Q = [0 0 1;1 0 0;0 1 0];
        vw = [-50 0];
    case 'armadillo'
        load('../ptclouds/armadillo.mat');
        load('../ptclouds/armadillo_patches.mat');
        Q = eye(3);
        vw = [-90 4];
    case 'children'
        load('../ptclouds/dancing_children.mat');
        load('../ptclouds/dancing_children_patches.mat');
        Q = [0 0 1;1 0 0;0 1 0];
        vw = [100 -1];
        % vw = [123 6];
    case 'pump'
        load('../ptclouds/pump_carter.mat');
        load('../ptclouds/pump_carter_patches.mat');
        Q = [0 0 1;1 0 0;0 1 0];
        vw = [78 3];
    case 'filigree'
        load('../ptclouds/filigree.mat');
%         load('../ptclouds/filigree_patches.mat');
        load('../ptclouds/filigree_patches_M35130.mat');
        Q = [0 0 1;1 0 0;0 1 0];
        vw = [90 0];
    case 'gargoyle'
        load('../ptclouds/gargoyle.mat');
        load('../ptclouds/gargoyle_patches.mat');
        Q = [0 0 1;-1 0 0;0 -1 0];
        vw = [-59 1];
    case 'bunny'
        load('../ptclouds/bunny_large.mat')
        load('../ptclouds/bunny_large_patches.mat')
%         load('../oldptclouds/bunny_h0.01.mat');
%         load('../oldptclouds/bunny_smooth_patches_M1447.mat');
        Q = [0 0 1;1 0 0;0 1 0];
%         Q = eye(3);
        vw = [93 4];
    case 'dragon'
        load('../ptclouds/stanford_dragon_fullres_face_nrmls.mat'); 
        [xc,ia,ic] = unique(xc,'rows'); nrml = nrml(ia,:);
%         load('../ptclouds/stanford_dragon_fullres_patches_N14400.mat');
        load('../ptclouds/stanford_dragon_patches_M32527.mat');
        Q = [0 0 1;1 0 0;0 1 0];
        vw = [60 15];
    case 'tori'
        load('../ptclouds/interlocked_tori_big.mat');
        load('../ptclouds/interlocked_tori_big_patches.mat');
        Q = eye(3);
        vw = [90 0];
    case 'hand'
        load('../ptclouds/laurent_hand.mat');
        load('../ptclouds/laurent_hand_patches.mat');
        Q = [0 0 1;1 0 0;0 1 0];
        vw = [90 0];
    case 'cantius_tooth'
        load('../ptclouds/cantius_tooth.mat');
        rng(6232020)
        [~,y] = kmeans(xc,200,'Replicates',5);
        Q = eye(3);
        vw = [-128 36];
    case 'mammoth_tooth'
        load('../ptclouds/mammoth_tooth.mat');
        rng(6232020)
        [~,y] = kmeans(xc,200,'Replicates',5);
        Q = eye(3);
        vw = [-128 36];
    case 'frog'
        load('../ptclouds/frog.mat');
        rng(6232020)
        [~,y] = kmeans(xc,100,'Replicates',5);
        Q = eye(3);
        vw = [0 90];
    otherwise
        load('../ptclouds/bunny_large.mat')
        load('../ptclouds/bunny_large_patches.mat')
        Q = [0 0 1;1 0 0;0 1 0];
        vw = [93 4];
end

%%
N = length(xc);
% tic
% [~,d] = knnsearch(xc,xc,'k',2);
% surface_area = N*mean(d(:,2))^2;
% y = util.PcCoarsen3D(xc,floor(N/15),surface_area);
% toc
tic
Nc = floor(N/40);
% Nc = 9000;
y = util.pcCoarsenPoissonDisk(xc,Nc);
% y = util.pcCoarsenWse(xc,Nc);
toc
[Nc length(y)]
x = xc;

tic
nrml = util.pcComputeNormals(x,10,2);
toc

%%
% id = randperm(N);
% id = id(1:floor(N/4));
% x = x(id,:);
% nrml = nrml(id,:);

%%

% Add noise to the normals
% id = 1:N;
% rng(1378271);
% nrml(id,:) = nrml(id,:) + 0.3*randn(size(nrml(id,:)));

% Make the normals unit vectors
nrml = nrml./(sqrt(sum(nrml.^2,2)));

% Rotate or shift the points
x = x*Q';
y = y*Q';
nrml = nrml*Q';
    
% Evaluation grid size
% mm = 56;
% mm = 512 - 7;
% mm = 436-7;
% mm = 384-7;
mm = 256-7;
% mm = 160;

% Regularization parameters:
reginfo.exactinterp = 1;
reginfo.nrmlreg = 0;
reginfo.nrmllambda = 1e-4;
reginfo.potreg = 0; 
reginfo.potlambda = 1e-2;

% reginfo.regularizationi = 2;  % Regularization for Homer used for MGM paper
% reginfo.lambdai = 1e1;

% Kernels to use
% RBFs to use
kernelinfo.phi = @(r) -r;
% kernelinfo.eta = @(r) r.^3;
% kernelinfo.zeta = @(r) 3*r;  kernelinfo.order = 2;
kernelinfo.eta = @(r) -r;
kernelinfo.zeta = @(r) -1./r; kernelinfo.order = 1;

plot3(x(:,1),x(:,2),x(:,3),'r.','MarkerSize',0.2)
hold on
plot3(y(:,1),y(:,2),y(:,3),'b.','MarkerSize',12);
hold off
daspect([1 1 1])

%%
% tic
% puminfo = cfpufit(x,nrml,y,kernelinfo,reginfo);
% toc
% 
% tic
% [potential,X,Y,Z] = cfpuval(puminfo,mm);
% toc

tic
[potential,X,Y,Z] = cfpurecon(x,nrml,y,kernelinfo,reginfo,mm);
toc

%%
clf
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
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
% view([90 15]);      % Buddha 
% view([-128 36]);    % cantius_tooth
% view([90 0]);     % Homer & Knot & Filigree
% view([100 -1]);   % Children
% view([123 6])     % Chidren zoom
% view([93 4]);     % Bunny
% view([78 3]);     % Carter
% view([-59 1]);      % Gargoyle
% view([-50 0]);      % Raptor
view(vw);
% camlight headlight
camlight('right','infinite')
% camlight('headlight','infinite');  % Bunny
% l = light('Position',[100 10 0],'Style','infinite');
% l2 = light('Position',[-100 0 0],'Style','infinite');
axis off
set(gcf,'Color',[1 1 1]*1)
% view([-40 20])
% view([-30 10]);
axis tight
ax = axis;


% axis([0.25 0.45 0 0.3 0.55 0.8]) % Gargoyle zoom
% axis([0.25 0.45 0 0.25 0.55 0.73]) % Gargoyle zoom
% axis([0.1 0.27 0.18 0.35 0.28 0.45]) % Children zoom

% %%
% fv = isosurface(X,Y,Z,potential,0);
% p = patch(fv);
% isonormals(X,Y,Z,potential,p);
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
% ww = nrml([a:,3);
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
%%
% verts = get(p, 'Vertices');
% faces = get(p, 'Faces');
% a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
% b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
% c = cross(a, b, 2);
% area = 1/2 * sum(sqrt(sum(c.^2, 2)));
% fprintf('\nThe surface area is %f\n\n', area);

%%

%%
% ns = zeros(length(y),1);
% for k=1:length(y)
%     ns(k) = puminfo.patchinfo(k).n;
% end
% histogram(ns)

% fittime = zeros(5,1);
% evaltime = fittime;
% 
% for j = 1:5
%     tic
%     puminfo = cfpufit(x,nrml,y,kernelinfo,reginfo);
%     fittime(j) = toc;
% 
%     tic
%     [potential,X,Y,Z] = cfpuval(puminfo,mm);
%     evaltime(j) = toc;
% end
% 
% fprintf('%s, fit time = %1.2f eval time = %1.2f\n\n',ptcloud,min(fittime),min(evaltime))