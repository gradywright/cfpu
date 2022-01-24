%% Set-up kernel struct used for the reconstruction

% Order 1 polyharmonic spline
kernelinfo.phi = @(r) -r;
kernelinfo.eta = @(r) -r;
kernelinfo.zeta = @(r) -1./r; kernelinfo.order = 1;

% Order 2 polyharmonic spline
% kernelinfo.phi = @(r) -r;
% kernelinfo.eta = @(r) r.^3;
% kernelinfo.zeta = @(r) 3*r;  kernelinfo.order = 2;

%% Homer
% Load the point cloud
load('../ptclouds/homer.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to exact interpolation and no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 256;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([90 0]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Trefoil Knot - no noise
% Load the point cloud
load('../ptclouds/trefoil_N23064.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to no regularization
reginfo.exactinterp = 1;
reginfo.nrmlreg = 0;
reginfo.potreg = 0; 
% Grid size
m = 150;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([90 0]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])


%% Trefoil Knot with noise - user specified smoothing of normals
% Load the point cloud
load('../ptclouds/trefoil_N23064.mat');
% Add noise to the normals
id = 1:length(ptcloud);
rng(1378271);
normals(id,:) = normals(id,:) + 0.3*randn(size(normals(id,:)));
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to exact interpolation and nrml regularization
reginfo.exactinterp = 1;
reginfo.nrmlreg = 1;        % User specified regularization
reginfo.nrmllambda = 1e-4;  % Ridge regression regularization parameter
reginfo.potreg = 0; 
% Grid size
m = 150;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);

% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([90 0]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Trefoil Knot with noise - GCV smoothing of normals
% Load the point cloud
load('../ptclouds/trefoil_N23064.mat');
% Add noise to the normals
id = 1:length(ptcloud);
rng(1378271);
% Make sure the normals unit vectors
normals(id,:) = normals(id,:) + 0.3*randn(size(normals(id,:)));
% Set regularization parameters to GCV regularization
reginfo.exactinterp = 1;
reginfo.nrmlreg = 2;
reginfo.potreg = 0; 
% Grid size
m = 150;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);

% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([90 0]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Interlocked Tori
% Load the point cloud
load('../ptclouds/interlocked_tori.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 128;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([90 15]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Simple Frog
% Load the point cloud
load('../ptclouds/frog.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 128;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([0 90]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Happy Buddha
% Load the point cloud
load('../ptclouds/happy_buddha.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 350;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([90 15]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Raptor head
% Load the point cloud
load('../ptclouds/raptor_head.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to exact interpolation and no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 256;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([-50 0]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Armadillo
% Load the point cloud
load('../ptclouds/armadillo.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to exact interpolation and no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 256;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([-90 4]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Pump carter
% Load the point cloud
load('../ptclouds/pump_carter.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to exact interpolation and no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 256;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([78 3]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Filigree
% Load the point cloud
load('../ptclouds/filigree.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to exact interpolation and no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 350;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([90 0]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Stanford Dragon
% Load the point cloud
load('../ptclouds/stanford_dragon.mat'); 
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to exact interpolation and no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 300;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([60 15]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Stanford Bunny
% Load the point cloud
load('../ptclouds/stanford_bunny.mat')
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to exact interpolation and regularization of
% normals (nrml) and the potential (pot)
reginfo.exactinterp = 1;
reginfo.nrmlreg = 1;        % User specified regularization of nrml
reginfo.nrmllambda = 1e-2;  % Ridge regression nrml regularization parameter
reginfo.potreg = 1;         % User specified regularization of pot
reginfo.potlambda = 1e-4;   % Ridge regression pot regularization parameter
% Grid size
m = 256;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);

% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([93 4]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Dancing children - no smoothing
% Load the point cloud
load('../ptclouds/dancing_children.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to exact interpolation and no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 300;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([100 -1]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Dancing children zoom
load('../ptclouds/dancing_children.mat');
vw = [123 6];

%% Gargoyle - no smoothing
% Load the point cloud
load('../ptclouds/gargoyle.mat');
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));
% Set regularization parameters to exact interpolation and no regularization
reginfo.exactinterp = 1; reginfo.nrmlreg = 0; reginfo.potreg = 0; 
% Grid size
m = 300;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([-59 1]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%%
% Make sure the normals unit vectors
normals = normals./(sqrt(sum(normals.^2,2)));

% Evaluation grid size
mm = 256;

% Regularization parameters:
reginfo.exactinterp = 1;
reginfo.nrmlreg = 0;
reginfo.nrmllambda = 1e-4;
reginfo.potreg = 0; 
reginfo.potlambda = 1e-2;

% plot3(x(:,1),x(:,2),x(:,3),'r.','MarkerSize',0.2)
% hold on
% plot3(y(:,1),y(:,2),y(:,3),'b.','MarkerSize',12);
% hold off
% daspect([1 1 1])

%%
% tic
% puminfo = cfpufit(x,nrml,y,kernelinfo,reginfo);
% toc
% 
% tic
% [potential,X,Y,Z] = cfpuval(puminfo,mm);
% toc

tic
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,kernelinfo,reginfo,mm);
toc

%%
clf
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view(vw);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1]*1)

% axis([0.25 0.45 0 0.3 0.55 0.8]) % Gargoyle zoom
% axis([0.25 0.45 0 0.25 0.55 0.73]) % Gargoyle zoom
% axis([0.1 0.27 0.18 0.35 0.28 0.45]) % Children zoom
