% clear all
% load('bunny_h0.01.mat')
% load('tooth8_h0.04.mat')
% load('frog');
load('happy_buddha_3.mat');
% load('homer.mat');

x = xc;
N = size(x,1);
% u = @(x) [2*x(:,1).*(4*x(:,1).^6-1), 2*x(:,2).*(4*x(:,2).^6-1), 2*x(:,3).*(4*x(:,3).^6-1)];
% nrml = u(x);
nrml = nrml./(sqrt(sum(nrml.^2,2)));

% Set up meshgrid for contour plot
[minx,maxx] = bounds(x);
mm = 60;
dx = (maxx-minx)/mm;
[X, Y, Z] = meshgrid(linspace(minx(1)-dx(1),maxx(1)+dx(1),mm),linspace(minx(2)-dx(2),maxx(2)+dx(2),mm),linspace(minx(3)-dx(3),maxx(3)+dx(3),mm));

% Add a point on the boundary of the shape 
% for setting the potential constant
xe = [X(:) Y(:) Z(:)];
m = length(xe);

A = 4*pi*max((sum((x-mean(x)).^2,2)));

n = 56;  % Number of nodes per patch
% patchRad = sqrt(A/4/pi*(n/N));   % Radius of the circular patches.
patchRad = sqrt(A/pi*(n/N));   % Radius of the circular patches.
q = 2;                           % Average number of patches a node belongs to
% M = ceil(A/(4*pi*patchRad^2));
M = ceil(A/(pi*patchRad^2));
patchRad = 2*patchRad;

%%
tic
rng(6232020)
[~,y] = kmeans(x,M,'Replicates',5);
toc

%%
plot3(x(:,1),x(:,2),x(:,3),'r.')
hold on
% plot3(y(:,1),y(:,2),y(:,3),'r.','MarkerSize',14)
for j=1:M
    [xx,yy,zz] = sphere(20);
    xx = xx*patchRad;
    yy = yy*patchRad;
    zz = zz*patchRad;
    xx = y(j,1)-xx;
    yy = y(j,2)-yy;
    zz = y(j,3)-zz;
    h = surf(xx,yy,zz);
    set(h,'FaceAlpha',0.0,'FaceColor',0.8*[0 0 1],'EdgeColor','k')
end
daspect([1 1 1])