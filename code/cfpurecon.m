function [potential,X,Y,Z] = cfpurecon(x,nrml,y,kernelinfo,reginfo,gridsize)
% CFPURECON Reconstructs a surface from a point cloud using the Curl-free
% Partition of Unity (CFPU) method.
%
% [P,X,Y,Z] = CFPURECON(X,NRML,Y,KERNELINFO,REGINFO,GRIDSIZE) reconstructs a 
% surface for the point cloud X using the normals NRML, and PU pathes Y.
% KERNELINFO contains the information about the curl-free kernel to use, and
% REGINFO specifies what regularization parameters should be used. GRIDSIZE
% specifies the size of the background grid to use for the isosurface (marching 
% cubes) extraction of the surface. This should be a 1-by-3 array [NX NY NZ]
% where NX, NY, NZ specifiy the size of the grid in the X, Y, and Z directions,
% respectively.  The larger these values, the crisper the surface will look, but
% the more time the code will take.
%
% A level surface can be obtained using the code
%       fv = isosurface(X,Y,Z,P,0);
%       ptch = patch(fv);
%       isonormals(X,Y,Z,P,ptch);
%       daspect([1 1 1])
%       view(3)
%
% The KERNELINFO structure should contain fields PHI, ETA, and ZETA. PHI is the
% scalar RBF used construct the curl-free kernel.  ETA should be 1/r (dPHI/dr)
% and ZETA = 1/r (dETA/dr).
%
% The REGINFO structure can contain following fields
% exactinterp: specfies whether the point cloud should exactly interpolate the
%              surface (if no regularization of the potential is used).
%          0 - no exact interpolation           
%          1 - enforce exact interpolation (default)
% nrmlreg: sets the regularization used in the fit of the normals
%          0 - no regularization (default)
%          1 - use ridge regression (smoothing splines) with a specified value
%              for the regularization parameter
%          2 - use ridge regression with the reg. parameter chosen using GCV
% nrmllambda: ridge regression parameter for the CF fit of the normals if
%             nrmlreg=1. Should be a value >= 0.
% potreg: sets the regularization used in the fit of the potential
%          0 - no regularization (default)
%          1 - use ridge regression (smoothing splines) with a specified value
%              for the regularization parameter
%          2 - use ridge regression with the reg. parameter chosen using GCV
% potlambda: ridge regression parameter for the fit of the potential to the
%            point cloud if potreg=1. Should be a value >= 0.
%
% Note that this code is faster than calling cfpufit and cfpuval, but the
% fitting data is not stored.
%
% see also CFPUFIT and CFPUVAL

% Copyright 2022 by Grady B. Wright

% Shift and scale the points to fit in [0,1]^3;
[minxx,maxxx] = bounds(x);
x = x - minxx;
x = x./max(maxxx-minxx);

% Shift and scale the centers
y = y - minxx;
y = y./max(maxxx-minxx);
M = size(y,1);

% Set up meshgrid for computing the potential
[minx,maxx] = bounds(x);

% Take each column of the patches to aid in better parfor performance.
patchx = y(:,1);
patchy = y(:,2);
patchz = y(:,3);

% Determine the patch radius
tree = createns(y);
[~,nn_dist] = knnsearch(tree,y,'k',2);
H = max(nn_dist(:,2));
delta = 1;
patchRad = (1 + delta)*H/2;

% Determine which nodes belong to which patch
[idx,nn_dist] = rangesearch(x,y,patchRad);
% idx = rangesearch(tree,x,patchRad);
% node_vec = cell(1,N);
% for j=1:N
%     id = idx{j};
%     for k = 1:length(id)
%         node_vec{id(k)} = [node_vec{id(k)} j];  
%     end
% end
% idx = node_vec;

% Regularization parameters
[exactinterp,nrmlreg,nrmllambda,nrmlschur,potreg,potlambda,trbl_id] = parseRegParams(reginfo);

% Radial kernels to use
eta = kernelinfo.eta;       % Curl-free 
zeta = kernelinfo.zeta;     % Curl-free 
phi = kernelinfo.phi;       % Exact interpolation of potential

% Polynomials to include
order = kernelinfo.order;

% Degree of the curl-free polynomial basis
if order == 1
    l = 3;
elseif order == 2
    l = 9;
else
    error('Curl-free polynomial degree not supported, choose 1 or 2');
end

% Set-up the background grid info
griddx = (maxx-minx)/gridsize;
griddx = max(griddx);
startx = (minx(1)-3*griddx);  endx = (maxx(1)+3*griddx);
starty = (minx(2)-3*griddx);  endy = (maxx(2)+3*griddx);
startz = (minx(3)-3*griddx);  endz = (maxx(3)+3*griddx);
xx = startx:griddx:endx;
yy = starty:griddx:endy;
zz = startz:griddx:endz;
[X, Y, Z] = meshgrid(xx,yy,zz);
[mmy,mmx,mmz] = size(X);
m = mmx*mmy*mmz;

% Variables for sparse storage
idxe_patch = cell(1,M);
patch_vec = cell(1,M);

% Matrices/vectors that are used over and over again
% zm = zeros(l);
% zv = zm(:,1);

% Initialize PU weight functions and CFPU potentials
Psi = cell(1,M);              % values for pum weight function on each patch
potential_local = cell(1,M);  % values for potential function

% % Array of structures for holding the patch data
% puminfo = struct('x',x,'y',y,'patchRad',patchRad,'idxinfo',[],'patchinfo',[],'kernelinfo',kernelinfo,'reginfo',reginfo);
% patchinfo = repmat(struct(...
%     'n',0,'coeffsx',[],'coeffsy',[],'coeffsz',[],'coeffsCorrection',[]),M,1);
% idxinfo = repmat(struct('id',[]),M,1);

opts.SYM = true;

% Loop over each patch and store local interpolant and Wendland function
parfor k = 1:M
    id = idx{k};
    h2 = max(nn_dist{k})^2;
    x_local = x(id,:);      % Grab nodes on patch
    xx_local = x_local(:,1).';
    xy_local = x_local(:,2).';
    xz_local = x_local(:,3).';
    n = size(x_local,1);        % Number of nodes on patch
    
    [CFP,P] = curlfreePoly(x_local,order);
    CFPt = CFP.';
    
    % Calculate the distance matrix between the nodes on each patch
    dx = (xx_local.' - xx_local);
    dy = (xy_local.' - xy_local);
    dz = (xz_local.' - xz_local);
    r = sqrt(dx.^2 + dy.^2 + dz.^2);

    % Pack the coefficient vector with the normal
    ui = nrml(idx{k},:);
    b = zeros(3*n+l,1);
    b(1:3:3*n,1) = ui(:,1);
    b(2:3:3*n,1) = ui(:,2);
    b(3:3:3*n,1) = ui(:,3);
    
    % Construct the interpolation matrix
    A = zeros(3*n+l);
    eta_temp = eta(r);
    zeta_temp = zeta(r);

    % Handle any potential singularities in zeta
    zeta_temp(1:n+1:end) = 0;
    
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
    
    A(1:3*n,3*n+1:end) = CFP;
    A(3*n+1:end,1:3*n) = CFPt;
    
    if ( nrmlreg ~= 2 )
        if ( nrmlreg == 1 )   % Manual regularization
            A(1:3*n,1:3*n) = A(1:3*n,1:3*n) + 3*n*nrmllambda*eye(3*n);
        elseif ( nrmlreg == 3 ) % Only regularize in local areas
            if any(trbl_id(idx{k}))
                A(1:3*n,1:3*n) = A(1:3*n,1:3*n) + 3*n*nrmllambda*eye(3*n);
            end
        end
        % Compute the coefficients
        if ( nrmlschur == 0 )  % Solve using standard Gaussian elimination
            coeffs = linsolve(A,b,opts);
            coeffsp = coeffs(3*n+1:end);
            coeffs = coeffs(1:3*n);
        else % Solve using Schur complement of the saddle point system
            A = A(1:3*n,1:3*n);
            b = b(1:3*n);
            coeffsp = pinv(CFPt*(A\CFP))*(CFPt*(A\b));
            coeffs = A\(b-CFP*coeffsp);
        end
    else % GCV regularization
        A = A(1:3*n,1:3*n);
        b = b(1:3*n);
        L = size(CFP,2);
        [F1,G1] = qr(CFP);
        F2 = F1(:,L+1:end);
        F1 = F1(:,1:L);
        G1 = G1(1:L,1:L);
        w1 = F1'*b;
        w2 = F2'*b;
        L = chol(F2'*A*F2);
        [U,D,~] = svd(L');
        D = diag(D);
        z = U'*w2;
        
        % Determine the parameter
        % [lam,~,flag,~] = fminbnd(@gcvCostFunction,-10,35,[],z,D,3*n);
        % lam = 3*n*exp(-lam);
        [lam,~,flag,~] = fminbnd(@gcvCostFunction,-10,35,[],z,D,3/h2);
        lam = 3/h2*exp(-lam);
        
        A = A + lam*eye(3*n);
        
        % Compute the coefficients
        coeffs = F2*(U*(z./(D.^2 + lam)));
        coeffsp = G1\(w1-F1'*(A*coeffs));
    end
    
    % Make things faster by setting temp variables.
    coeffsx = coeffs(1:3:3*n).';
    coeffsy = coeffs(2:3:3*n).';
    coeffsz = coeffs(3:3:3*n).';    

    temp_potential_nodes = sum(eta_temp.*(dx.*coeffsx + dy.*coeffsy + dz.*coeffsz),2) + P*coeffsp;
    
    % Use a scalar RBF fit of the residual to correct the potential
    if ( exactinterp ) 
        % Append degree 0 polynomial
        P = ones(n,1);
        A = ones(n+1);
        A(1:n,1:n) = phi(r);
        A(end) = 0;
        b = [temp_potential_nodes;0];
        if ( potreg ~= 2 )
            if ( potreg == 1 )
                A(1:n,1:n) = A(1:n,1:n) + n*potlambda*eye(n);
            elseif ( potreg == 3 )
                if any(trbl_id(idx{k}))
                    A(1:n,1:n) = A(1:n,1:n) + n*potlambda*eye(n);
                end
            end
            coeffs_correction = linsolve(A,b,opts);
        else
            L = size(P,2);
            b = b(1:n);
            A = A(1:n,1:n);
            [F1,G1] = qr(P);
            F2 = F1(:,L+1:end);
            F1 = F1(:,1:L);
            G1 = G1(1:L,1:L);
            w1 = F1'*b;
            w2 = F2'*b;
            L = chol(F2'*A*F2);
            [U,D,~] = svd(L');
            D = diag(D);
            z = U'*w2;
            
            % Determine the parameter
            % [lam,~,flag,~] = fminbnd(@gcvCostFunction,-10,35,[],z,D,n);
            % lam = n*exp(-lam);
            [lam,~,flag,~] = fminbnd(@gcvCostFunction,-10,35,[],z,D,1/h2);
            lam = (1/h2)*exp(-lam);
            A = A + lam*eye(n);

            % Compute the coefficients
            temp = F2*(U*(z./(D.^2 + lam)));
            coeffs_correction = [temp;G1\(w1-F1'*(A*temp))];
        end
    else  % Use a linear fit of the residual to correct the potential
        P = [P(:,1:3) ones(n,1)];
        coeffs_correction = P\temp_potential_nodes;
    end    
    
    coeffs_correction_const = coeffs_correction(end);
    coeffs_correction = coeffs_correction(1:end-1);
    
    % Determine the local evaluation points using some fancy indexing for grids    
    ix = round((patchx(k)-startx)/griddx)+1;
    iy = round((patchy(k)-starty)/griddx)+1;
    iz = round((patchz(k)-startz)/griddx)+1;
    factor = round(patchRad/griddx);
    ixs = max(ix-factor,1):min(ix+factor,mmx);
    iys = max(iy-factor,1):min(iy+factor,mmy);
    izs = max(iz-factor,1):min(iz+factor,mmz);
    xx = (startx+(ixs-1)*griddx);
    yy = (starty+(iys-1)*griddx).';
    zz = shiftdim(startz+(izs-1)*griddx,-1);
    De = (patchx(k)-xx).^2 + (patchy(k)-yy).^2 + (patchz(k)-zz).^2;
    id = De(:) < patchRad^2;
    ixs2 = repmat(ixs,[length(yy) 1 length(zz)]);
    iys2 = repmat(iys.',[1 length(xx) length(zz)]);
    izs2 = repmat(shiftdim(izs,-1),[length(yy) length(xx) 1]);
    temp_idg = (iys2 + (ixs2 - 1).*mmy) + (izs2-1)*(mmx*mmy);
    temp_idg = temp_idg(id);    
    De = sqrt(De(id));
    idxe_patch{k} = temp_idg;

    % Calculate the weight function on each patch center and store it
    Psi{k} = weight(De,patchRad,0);

    % Local grid points in the patch
    mlocalx = length(xx);
    mlocaly = length(yy);
    mlocalz = length(zz);
    xx = repmat(xx,[mlocaly 1 mlocalz]);
    yy = repmat(yy,[1 mlocalx mlocalz]);
    zz = repmat(zz,[mlocaly mlocalx 1]);
    xe_local = [xx(id) yy(id) zz(id)];
    
    mm = size(xe_local,1);         % Number of evaluation points
    if mm == 0
        continue
    end

    % Batch up the evaluation points to make things faster.
    batch_sz = ceil(100^2/n);
    temp_potential = zeros(mm,1);
    potential_correction = zeros(mm,1);
    for j = 1:batch_sz:mm
        idb = j:min((j+batch_sz-1),mm);
        xe_local_batch = xe_local(idb,:);
        dx = (xe_local_batch(:,1) - xx_local);
        dy = (xe_local_batch(:,2) - xy_local);
        dz = (xe_local_batch(:,3) - xz_local);
        r = sqrt(dx.^2 + dy.^2 + dz.^2);
        
        [~,P] = curlfreePoly(xe_local_batch,order);
        
        temp_potential(idb) = sum(eta(r).*(dx.*coeffsx + dy.*coeffsy + dz.*coeffsz),2) + P*coeffsp;

        % Potential
        if exactinterp
            potential_correction(idb) = phi(r)*coeffs_correction +  coeffs_correction_const;
        else
            potential_correction(idb) = P(:,1:3)*coeffs_correction + coeffs_correction_const;
        end
    end
    
    % Correct the potential function
    potential_local{k} = temp_potential - potential_correction;    
    patch_vec{k} = k*ones(mm,1);
end

% Evaluate the potentials over each patch and combine with weight functions
patch_vec = vertcat(patch_vec{:});
idxe_vec = vertcat(idxe_patch{:});
Psi_sum = sum(sparse(idxe_vec,patch_vec,vertcat(Psi{:}),m,M),2);
for k = 1:M
    potential_local{k} = potential_local{k}.*(Psi{k}./Psi_sum(idxe_patch{k}));
end
% Compute the potential function
temp = sum(sparse(idxe_vec,patch_vec,vertcat(potential_local{:}),m,M),2);
% Put nans where the evaluation points are not included in any patches
[i,~] = find(Psi_sum);
potential = nan(m,1);
potential(i) = temp(i);
potential = reshape(potential,size(X));

end

function [exactinterp,nrmlreg,nrmllam,nrmlschur,potreg,potlam,trbl_id] = parseRegParams(reginfo)

% Default values for the regularization parameters
exactinterp = 1;
nrmlreg = 0;
nrmllam = 0;
nrmlschur = 0;
trbl_id = cell(0);
potreg = 0;
potlam = 0;

if ( isfield(reginfo,'exactinterp') )
    exactinterp = reginfo.exactinterp;
end

if ( isfield(reginfo,'nrmlreg') )
    nrmlreg = reginfo.nrmlreg;
end

if ( isfield(reginfo,'nrmlschur') )
    nrmlschur = reginfo.nrmlschur;
end

if ( isfield(reginfo,'nrmllambda') )
    nrmllam = reginfo.nrmllambda;
end

if ( nrmlreg == 3 )
    if ( isfield(reginfo,'trbl_id') )
        trbl_id = reginfo.trbl_id;
    end
end

if ( isfield(reginfo,'potreg') )
    potreg = reginfo.potreg;
end

if ( isfield(reginfo,'potlambda') )
    potlam = reginfo.potlambda;
end

end