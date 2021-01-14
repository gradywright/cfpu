% load('../ptclouds/bunny_large.mat')
% load('../ptclouds/tooth8_h0.10.mat')
% load('../ptclouds/frog.mat');
% load('../ptclouds/happy_buddha.mat');
% load('../ptclouds/stanford_dragon_fullres_face_nrmls.mat'); [xc,ia,ic] = unique(xc,'rows'); nrml = nrml(ia,:);
% load('../ptclouds/happy_buddha_highres.mat');
% load('../ptclouds/homer.mat');
% load('../ptclouds/armadillo.mat');
% load('../ptclouds/mammoth_tooth.mat');
% load('../ptclouds/cantius_tooth.mat');
% load('../ptclouds/interlocked_tori_big.mat');
% load('../ptclouds/laurent_hand.mat');
% load('../ptclouds/raptor.mat');
load('../ptclouds/raptor_head.mat');
% load('../ptclouds/pump_carter.mat');
% load('../ptclouds/filigree.mat');
% load('../ptclouds/gargoyle.mat');
% load('../ptclouds/dancing_children.mat');
% load('../ptclouds/trefoil_N23064.mat');
% load('../ptclouds/trefoil_N6144.mat');

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
% nrml(id,:) = nrml(id,:) + 0.3*randn(size(nrml(id,:)));

nrml = nrml./(sqrt(sum(nrml.^2,2)));

Q = [0 0 1;1 0 0;0 1 0];
% Q = [0 0 1;-1 0 0;0 -1 0];
% Q = eye(3);
x = x*Q';
nrml = nrml*Q';

% Set up meshgrid for contour plot
[minxx,maxxx] = bounds(x);

x = x - minxx;
x = x./max(maxxx-minxx);

[minx,maxx] = bounds(x);

% mm = 56;
% mm = 512 - 7;
mm = 436-7;
% mm = 384-7;
% mm = 256-7;
% mm = 160;
dx = (maxx-minx)/mm;
% [X, Y, Z] = meshgrid(linspace(minx(1)-dx(1),maxx(1)+dx(1),mm),linspace(minx(2)-dx(2),maxx(2)+dx(2),mm),linspace(minx(3)-dx(3),maxx(3)+dx(3),mm));
% [X, Y, Z] = meshgrid(linspace(minx(1)-4*dx(1),maxx(1)+4*dx(1),mm),linspace(minx(2)-4*dx(2),maxx(2)+4*dx(2),mm),linspace(minx(3)-4*dx(3),maxx(3)+4*dx(3),mm));
dx = max(dx);
xx = (minx(1)-3*dx):dx:(maxx(1)+3*dx);
yy = (minx(2)-3*dx):dx:(maxx(2)+3*dx);
zz = (minx(3)-3*dx):dx:(maxx(3)+3*dx);

% minx2 = [0.2 0 0.4];    % Gargoyle zoom
% maxx2 = [0.5 0.4 0.9];
% minx2 = [0.1 0.15 0.25];    % Children zoom
% maxx2 = [0.3 0.35 0.45];
% mm2 = 256-7;
% dx2 = (maxx2-minx2)/mm2;
% xx = (minx2(1)-3*dx2):dx2:(maxx2(1)+3*dx2);
% yy = (minx2(2)-3*dx2):dx2:(maxx2(2)+3*dx2);
% zz = (minx2(3)-3*dx2):dx2:(maxx2(3)+3*dx2);

[X, Y, Z] = meshgrid(xx,yy,zz);
xe = [X(:) Y(:) Z(:)];
m = length(xe);

% % Bounding box for troublsome points for dancing children
% trbl_bbox_min = [0.15 0.22 0.31];
% trbl_bbox_max = [0.29 0.28 0.36];
% trbl_id = all((x >= trbl_bbox_min) & (x <= trbl_bbox_max),2);

% % Bounding box for troublsome points for Gargoyle
% trbl_bbox_min = [0.36 0 0.63];
% trbl_bbox_max = [0.4 0.055 0.65];
% trbl_id = all((x >= trbl_bbox_min) & (x <= trbl_bbox_max),2);

% Number of nodes per patch
% n = 100;  % Homer
% n = 100;
% n = 200;
% n = 50;
% A = 4*pi*max((sum((x-mean(x)).^2,2)));
% patchRad = sqrt(A/pi*(n/N));   % Radius of the circular patches.
% q = 2;                         % Average number of patches a node belongs to
% M = ceil(A/(pi*patchRad^2))
% patchRad = 1.25*patchRad
%%
usekMeans = 0;

if usekMeans
    tic
    rng(6232020)
    [~,y] = kmeans(x,M,'Replicates',5);
    toc
else
%     load('../ptclouds/happy_buddha_highres_patches.mat');
%     load('../ptclouds/bunny_large_patches.mat');
%     load('../ptclouds/stanford_dragon_fullres_patches_N14400.mat');
%     load('../ptclouds/armadillo_patches.mat');
%     load('../ptclouds/happy_buddha_highres_patches.mat');
%     load('../ptclouds/interlocked_tori_big_patches.mat');
%     load('../ptclouds/laurent_hand_patches.mat');
%     load('../ptclouds/gargoyle_patches.mat');
%     load('../ptclouds/dancing_children_patches.mat');
%     load('../ptclouds/trefoil_patches_M864.mat');
%     load('../ptclouds/homer_patches.mat');
    load('../ptclouds/raptor_head_patches.mat');
%     load('../ptclouds/raptor_patches.mat');
%     load('../ptclouds/pump_carter_patches.mat');
%     load('../ptclouds/filigree_patches.mat')
    y = y*Q';
    
    y = y - minxx;
    y = y./max(maxxx-minxx);

    M = length(y);
end
tree = createns(y);
% patchRad = 2.5/sqrt(M);
[~,nn_dist] = knnsearch(tree,y,'k',2);
H = max(nn_dist(:,2));
% delta = 2/3;
delta = 1;
patchRad = (1 + delta)*H/2;

%%
% Determine which nodes belong to which patch
tic
idx = rangesearch(x,y,patchRad);
% idx = rangesearch(tree,x,patchRad);
% node_vec = cell(1,N);
% for j=1:N
%     id = idx{j};
%     for k = 1:length(id)
%         node_vec{id(k)} = [node_vec{id(k)} j];  
%     end
% end
% idx = node_vec;
toc

% Determine which evaluation nodes belong to which patch
tic
% [idxe,De] = rangesearch(xe,y,patchRad);
[eval_vec,temp_De] = rangesearch(tree,xe,patchRad);
idxe = cell(M,1);
% De = cell(M,1);
% temp_length = cellfun(@length,eval_vec);
% id_length = find(id > 0);
for j=1:m
    id = eval_vec{j};
    for k = 1:length(id)
        idxe{id(k)} = [idxe{id(k)} j];
%         De{id(k)} = [De{id(k)} temp_De{j}];
    end
end
toc

%%
% Fitting parameters:
lambda = 1e-3;
regularization = 0;
schurcmplmnt = 0;
exactinterp = 1;
regularizationi = 0;
lambdai = 1e-4;

% RBFs to use
phi = @(r) -r;
% eta = @(r) r.^3;
% zeta = @(r) 3*r;  order = 2;
eta = @(r) -r;
zeta = @(r) -1./r; order = 1;

% Variables for sparse storage
idxe_vec = [idxe{:}]';
patch_vec = cell(1,M);

% Perform RBF locally for each patch and store
su_local = cell(1,M); % matrix for local interpolants on each patch
sv_local = cell(1,M); % matrix for local interpolants on each patch
Psi = cell(1,M);     % matrix for Wendland function on each patch
Psi_grad = cell(1,M);     % matrix for Wendland function on each patch
potential_local = cell(1,M);     % matrix for potential function

lams = zeros(M,1);
flags = lams;
ns = zeros(M,1);
% Loop over each patch and store local interpolant and Wendland function
parfor k = 1:M
%     warning('off','MATLAB:nearlySingularMatrix')
    if mod(k,40) == 0
        fprintf('k = %d out of M = %d\n',k,M)
    end
    x_local = x(idx{k},:); % Grab nodes on patch
    n = size(x_local,1);    % Number of nodes on patch
    ns(k) = n;
    xe_local = xe(idxe{k},:);  % Evaluation points on each patch
    mm = size(xe_local,1);     % Number of evaluation points
    if mm == 0
        continue
    end
    
    [CFP,P] = curlfreePoly(x_local,order);

    % Least squares
    if regularization == 3
%         nlsq = min(max(ceil(n/4),10),n);
        nlsq = min(10,n);
        [~,xc_local] = kmeans(x_local,nlsq,'Replicates',5);
        [CFPt,~] = curlfreePoly(xc_local,order);
        CFPt = CFPt.';
    else
        xc_local = x_local;
        nlsq = n;
        CFPt = CFP.';
    end
    
    % Calculate the distance matrix between the nodes on each patch
    dx = (x_local(:,1) - xc_local(:,1).');
    dy = (x_local(:,2) - xc_local(:,2).');
    dz = (x_local(:,3) - xc_local(:,3).');
    r = sqrt(dx.^2 + dy.^2 + dz.^2);

    %Vector field coefficient vector
    ui = nrml(idx{k},:);
    b = zeros(3*n,1);
    b(1:3:3*n,1) = ui(:,1);
    b(2:3:3*n,1) = ui(:,2);
    b(3:3:3*n,1) = ui(:,3);
    
    %Construct A
    A = zeros(3*n,3*nlsq);
    eta_temp = eta(r);
    zeta_temp = zeta(r);
    if regularization ~= 3
        zeta_temp(1:n+1:end) = 0;
    else
        zeta_temp(r == 0) = 0;
    end
    dphi_xx = zeta_temp.*dx.^2+eta_temp;
    dphi_yy = zeta_temp.*dy.^2+eta_temp;
    dphi_zz = zeta_temp.*dz.^2+eta_temp;
    
    dphi_xy = zeta_temp.*dx.*dy;
    dphi_xz = zeta_temp.*dx.*dz;
    dphi_yz = zeta_temp.*dy.*dz;
        
    A(1:3:3*n,1:3:3*nlsq) = dphi_xx;
    A(1:3:3*n,2:3:3*nlsq) = dphi_xy;
    A(1:3:3*n,3:3:3*nlsq) = dphi_xz;
    
    A(2:3:3*n,1:3:3*nlsq) = dphi_xy;
    A(2:3:3*n,2:3:3*nlsq) = dphi_yy;
    A(2:3:3*n,3:3:3*nlsq) = dphi_yz;
    
    A(3:3:3*n,1:3:3*nlsq) = dphi_xz;
    A(3:3:3*n,2:3:3*nlsq) = dphi_yz;
    A(3:3:3*n,3:3:3*nlsq) = dphi_zz;
    
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
        L = size(CFP,2);
        [F1,G1] = qr(CFP);
        F2 = F1(:,L+1:end);
        F1 = F1(:,1:L);
        G1 = G1(1:L,1:L);
        w1 = F1'*b;
        w2 = F2'*b;
        L = chol(F2'*A*F2);
        [U,D,V] = svd(L');
        D = diag(D);
        z = U'*w2;
        
        % Determine the parameter
        [lam,fval,flag,output] = fminbnd(@gcv_cost_phs3,-10,35,[],z,D,3*n);
        lams(k) = 3*n*exp(-lam);
        flags(k) = flag;
        A = A + lams(k)*eye(3*n);
        
        % Compute the coefficients
        coeffs = F2*(U*(z./(D.^2 + lams(k))));
        coeffsp = G1\(w1-F1'*(A*coeffs));       
    else
        if regularization == 2   % Manual regularization
            A = A + 3*n*lambda*eye(3*n);
        elseif regularization == 4
            if any(trbl_id(idx{k}))
                A = A + 3*n*lambda*eye(3*n);
            end
        end
        % Compute the coefficients
        l = size(CFP,2);
        if schurcmplmnt == 0
            coeffs = [A CFP;CFPt zeros(l)]\[b;zeros(l,1)];
%             coeffs = pinv([A CFP;CFPt zeros(l)])*[b;zeros(l,1)];
            coeffsp = coeffs(3*nlsq+1:end);
            coeffs = coeffs(1:3*nlsq);        
        else
            coeffsp = pinv(CFPt*(A\CFP))*(CFPt*(A\b));
            coeffs = A\(b-CFP*coeffsp);
        end
    end
            
    temp_potential_nodes = sum(eta_temp.*(dx.*coeffs(1:3:3*nlsq)' + dy.*coeffs(2:3:3*nlsq)' + dz.*coeffs(3:3:3*nlsq)'),2) + P*coeffsp;
    
    if exactinterp
        P = ones(n,1);
        dx = (x_local(:,1) - x_local(:,1).');
        dy = (x_local(:,2) - x_local(:,2).');
        dz = (x_local(:,3) - x_local(:,3).');
        r = sqrt(dx.^2 + dy.^2 + dz.^2);
        A = phi(r);
        b = temp_potential_nodes;
        if regularizationi == 1
            L = size(P,2);
            [F1,G1] = qr(P);
            F2 = F1(:,L+1:end);
            F1 = F1(:,1:L);
            G1 = G1(1:L,1:L);
            w1 = F1'*b;
            w2 = F2'*b;
            L = chol(F2'*A*F2);
            [U,D,V] = svd(L');
            D = diag(D);
            z = U'*w2;

            % Determine the parameter
            [lam,fval,flag,output] = fminbnd(@gcv_cost_phs3,-10,35,[],z,D,n);
            lams(k) = n*exp(-lam);
            flags(k) = flag;
            A = A + lams(k)*eye(n);
            % Compute the coefficients
            temp = F2*(U*(z./(D.^2 + lams(k))));
            coeffs_correction = [temp;G1\(w1-F1'*(A*temp))];
        else
            if regularizationi == 2
                A = A + n*lambdai*eye(n);
            elseif regularizationi == 4
                if any(trbl_id(idx{k}))
                    A = A + n*lambdai*eye(n);
                end
            end
            coeffs_correction = [A P;P.' 0]\[temp_potential_nodes;0];
        end
    else
        P = [ones(n,1) P(:,1:3)];
%         P = ones(n,1);
        coeffs_correction = P\temp_potential_nodes;
    end
        
    % Generate the distance matrix between the patch nodes and global evaluation points    
    batch_sz = ceil(800^2/n);
    temp_potential = zeros(mm,1);
    potential_correction = zeros(mm,1);
    for j = 1:batch_sz:mm
        idb = j:min((j+batch_sz-1),mm);
        mmb = length(idb);
        dx = (xe_local(idb,1) - xc_local(:,1)');
        dy = (xe_local(idb,2) - xc_local(:,2)');
        dz = (xe_local(idb,3) - xc_local(:,3)');
        r = sqrt(dx.^2 + dy.^2 + dz.^2);
        
        [~,P] = curlfreePoly(xe_local(idb,:),order);
        
        temp_potential(idb) = sum(eta(r).*(dx.*coeffs(1:3:3*nlsq)' + dy.*coeffs(2:3:3*nlsq)' + dz.*coeffs(3:3:3*nlsq)'),2) + P*coeffsp;

        % Potential
        if exactinterp
            dx = (xe_local(idb,1) - x_local(:,1)');
            dy = (xe_local(idb,2) - x_local(:,2)');
            dz = (xe_local(idb,3) - x_local(:,3)');
            r = sqrt(dx.^2 + dy.^2 + dz.^2);
            potential_correction(idb) = [phi(r) ones(mmb,1)]*coeffs_correction;
        else
            potential_correction(idb) = [ones(mmb,1) P(:,1:3)]*coeffs_correction;
%             potential_correction(idb) = ones(mmb,1)*coeffs_correction;
        end
    end

    % Calculate the weight function on each patch center and store it
%     temp = De{k}.'/patchRad;
    De = sqrt(sum((y(k,:)-xe_local).^2,2));
%     Psi{k} = (max(1-De/patchRad,0));
%     Psi{k} = (max(1-temp,0).^4).*(4*temp+1);
%     Psi{k} = weight(De{k}.',patchRad,0);
    Psi{k} = weight(De,patchRad,0);
    
    % Correct the potential function
    potential_local{k} = temp_potential - potential_correction;
    
    patch_vec{k} = k*ones(mm,1);
end

%
% This code can be slow and require too much memory. The specific
% issue is with the line w = Psi_eval./Psi_sum;
%
% tic
% % Create sparse matrices for local interpolants and Wendland functions
% Psi_eval = sparse([idxe{:}]',patch_vec,vertcat(Psi{:}),m,M);
% 
% % Create weight matrix
% Psi_sum=sum(Psi_eval,2);
% w = Psi_eval./Psi_sum;
% % Sparse matrices for potential function and Psi_Qgrad
% potential_loc_eval = sparse([idxe{:}]',patch_vec,vertcat(potential_local{:}),m,M);
% g = full(sum((potential_loc_eval).*w,2));
% g = reshape(g,size(X));
% toc

%
% This code does the same thing as above but in a much more memory
% efficient manner.
patch_vec = vertcat(patch_vec{:});
Psi_sum = sum(sparse(idxe_vec,patch_vec,vertcat(Psi{:}),m,M),2);
for k = 1:M
    potential_local{k} = potential_local{k}.*(Psi{k}./Psi_sum(idxe{k}));
end
% Compute the potential function
temp = sum(sparse(idxe_vec,patch_vec,vertcat(potential_local{:}),m,M),2);
% Put nans where the evaluation points are not included in any patches
[i,j] = find(Psi_sum);
g = nan(m,1);
g(i) = temp(i);
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
% view([90 15]);   % Buddha 
% view([-128 36]);    % cantius_tooth
% view([90 0]);     % Knot & Filigree
% view([100 -1]);   % Children
% view([123 6])     % Chidren zoom
% view([93 4]);     % Bunny
% view([78 3]);     % Carter
% view([-59 1]);      % Gargoyle
view([-50 0]);      % Raptor
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

