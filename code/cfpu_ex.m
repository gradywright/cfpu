%% Homer
% Load the point cloud
load('../ptclouds/homer.mat');
% Make sure the normals unit vectors
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 256;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
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
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 150;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
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
normals = normals./sqrt(sum(normals.^2,2));
% Set the kernel to order 1 polyharmonic spline (the default)
kernel.phi = @(r) -r; kernel.eta = @(r) -r; kernel.zeta = @(r) -1./r; kernel.order = 1;
% Set regularization parameters to exact interpolation and nrml regularization
regularization.exactinterp = 1;
regularization.nrmlreg = 1;        % User specified regularization
regularization.nrmllambda = 1e-4;  % Ridge regression regularization parameter
regularization.potreg = 0; 
% Grid size
m = 150;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m,kernel,regularization);

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
regularization.exactinterp = 1;
regularization.nrmlreg = 2;
regularization.potreg = 0; 
% Grid size
m = 150;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m,kernel,regularization);

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
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 128;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
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
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 128;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
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
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 350;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
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
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 256;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
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
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 256;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
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
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 256;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
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
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 350;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
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
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 300;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
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
normals = normals./sqrt(sum(normals.^2,2));
% Set the kernel to order 1 polyharmonic spline (the default)
kernel.phi = @(r) -r; kernel.eta = @(r) -r; kernel.zeta = @(r) -1./r; kernel.order = 1;
% Set regularization parameters to exact interpolation and regularization of
% normals (nrml) and the potential (pot)
regularization.exactinterp = 1;
regularization.nrmlreg = 1;        % User specified regularization of nrml
regularization.nrmllambda = 1e-2;  % Ridge regression nrml regularization parameter
regularization.potreg = 1;         % User specified regularization of pot
regularization.potlambda = 1e-4;   % Ridge regression pot regularization parameter
% Grid size
m = 256;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m,kernel,regularization);

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
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 300;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([100 -1]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Gargoyle - no smoothing
% Load the point cloud
load('../ptclouds/gargoyle.mat');
% Make sure the normals unit vectors
normals = normals./sqrt(sum(normals.^2,2));
% Grid size
m = 300;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([-59 1]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

%% Example of reading a point cloud from a delimited text file
% Read in the points
ptcloud = load('../ptclouds/bozbezbozzel.txt');
% Remove any duplicate points from the point cloud
[ptcloud,ia,ic] = unique(ptcloud,'rows');
% Generate the normals using VCGLIb (use 10 neighbors and no smoothing)
normals = util.pcComputeNormals(ptcloud,10,0);
% Generate the patches using VCGLIb (use N/30 as rough estimate for the number)
patches = util.pcCoarsenPoissonDisk(ptcloud,ceil(length(ptcloud)/10));
% Grid size
m = 300;
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,m);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0);
p = patch(fv);
isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([62 0]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1])

