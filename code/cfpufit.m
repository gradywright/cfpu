function puinfo = cfpufit(x,nrml,y,kernelinfo,reginfo)
% CFPUFIT Fits a point cloud using the Curl-free Partition of Unity (CFPU)
% method.
%
% PUINFO = CFPUFIT(X,NRML,Y) Determines the expansion
% coefficients (i.e. the fit) of the surface that fits the point cloud X using
% the normals NRML, and PU pathes Y.  
%
% PUINFO = CFPUFIT(X,NRML,Y,KERNEL) Specifies what kernel should be used for
% constructing the curl-free RBF and the scalar correction field. The KERNEL
% structure should contain fields PHI, ETA, and ZETA. PHI is the scalar RBF used
% construct the curl-free kernel. ETA should be 1/r (dPHI/dr) and ZETA = 1/r
% (dETA/dr). The default is a order 1 polyharmonic spline for the curl-free RBF
% and an order 1 polyharmonic spline for the scalar correction. This is given by
% the code:
%       kernel.phi = @(r) -r;
%       kernel.eta = @(r) -r;
%       kernel.zeta = @(r) -1./r; 
%       kernel.order = 1;
% The code for using an order 2 polyharmonic spline for the curl-free RBF and
% order 2 polyharmonic spline for the scalar correction is given by the code:
%       kernel.phi = @(r) r.^3;
%       kernel.eta = @(r) r.^3;
%       kernel.zeta = @(r) 3*r;  
%       kernel.order = 2;
%
% Note that this code should be used when the same fit of the point cloud is
% to be reused.  Otherwise, one should just call cfpurecon directly
%
% see also CFPUVAL and CFPURECON

% Copyright 2022 by Grady B. Wright

% Check input arguments
if ( nargin < 3 )
    error('Not enough input arguments.');
elseif ( nargin == 3 )
    % Default kernel is order 1 polyharmonic spline
    kernelinfo.phi = @(r) -r;
    kernelinfo.eta = @(r) -r;
    kernelinfo.zeta = @(r) -1./r; 
    kernelinfo.order = 1;
    % Default is no regularization
    reginfo.exactinterp = 1;
elseif ( nargin == 4 )
    % Default is no regularization
    reginfo.exactinterp = 1;
end

% Shift and scale the points to fit in [0,1]^3;
[minxx,maxxx] = bounds(x);
x = x - minxx;
x = x./max(maxxx-minxx);

% Shift and scale the centers
y = y - minxx;
y = y./max(maxxx-minxx);
M = size(y,1);

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

% Matrices/vectors that are used over and over again
% zm = zeros(l);
% zv = zm(:,1);

% Array of structures for holding the patch data
puinfo = struct('x',x,'y',y,'patchRad',patchRad,'idxinfo',[],'patchinfo',[],'kernelinfo',kernelinfo,'reginfo',reginfo);
patchinfo = repmat(struct(...
    'n',0,'coeffsx',[],'coeffsy',[],'coeffsz',[],'coeffsCorrection',[]),M,1);
idxinfo = repmat(struct('id',[]),M,1);

opts.SYM = true;

% Loop over each patch and store local interpolant/approximant coefficients and
% PU weight function
parfor k = 1:M
    id = idx{k};
    h2 = max(nn_dist{k})^2;
    idxinfo(k).id = id;
    x_local = x(id,:);          % Grab nodes on patch
    xx_local = x_local(:,1).';
    xy_local = x_local(:,2).';
    xz_local = x_local(:,3).';
    
    n = size(x_local,1);        % Number of nodes on patch
    patchinfo(k).idx = idx{k};
    patchinfo(k).n = n;
    
    [CFP,P] = util.curlfreePoly(x_local,order);  % Calculate the CF polynomials for the patch nodes
    CFPt = CFP.';
    
    % Calculate the distance matrix between the nodes on each patch
    dx = (xx_local.' - xx_local);
    dy = (xy_local.' - xy_local);
    dz = (xz_local.' - xz_local);
    r = sqrt(dx.^2 + dy.^2 + dz.^2);
    
    % Pack the right hand side vector with the normals
    ui = nrml(idx{k},:);
    b = zeros(3*n+l,1);
    b(1:3:3*n,1) = ui(:,1);
    b(2:3:3*n,1) = ui(:,2);
    b(3:3:3*n,1) = ui(:,3);
    
    % Construct the interpolation matrix
    A = zeros(3*n+l);
    eta_temp = eta(r);
    zeta_temp = zeta(r);

    % Handle any potential (removable) singularities in zeta
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
%             A(1:3*n,1:3*n) = A(1:3*n,1:3*n) + 3*n*nrmllambda*eye(3*n);
              A(1:3*n,1:3*n) = A(1:3*n,1:3*n) + 3/h2*nrmllambda*eye(3*n);
        elseif ( nrmlreg == 3 ) % Only regularize in local areas
            if any(trbl_id(idx{k}))
%                 A(1:3*n,1:3*n) = A(1:3*n,1:3*n) + 3*n*nrmllambda*eye(3*n);
                A(1:3*n,1:3*n) = A(1:3*n,1:3*n) + 3/h2*nrmllambda*eye(3*n);
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
        % [lam,~,flag,~] = fminbnd(@util.gcvCostFunction,-10,35,[],z,D,3*n);
        % lam = 3*n*exp(-lam);
        [lam,~,flag,~] = fminbnd(@util.gcvCostFunction,-10,35,[],z,D,3/h2);
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
    
    patchinfo(k).coeffsp = coeffsp;
    patchinfo(k).coeffsx = coeffsx;
    patchinfo(k).coeffsy = coeffsy;
    patchinfo(k).coeffsz = coeffsz;
    
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
                % A(1:n,1:n) = A(1:n,1:n) + n*potlambda*eye(n);
                A(1:n,1:n) = A(1:n,1:n) + (1/h2)*potlambda*eye(n);
            elseif ( potreg == 3 )
                if any(trbl_id(idx{k}))
                    % A(1:n,1:n) = A(1:n,1:n) + n*potlambda*eye(n);
                    A(1:n,1:n) = A(1:n,1:n) + (1/h2)*potlambda*eye(n);
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
    patchinfo(k).coeffs_correction = coeffs_correction;
end

puinfo.patchinfo = patchinfo;
puinfo.idxinfo = idxinfo;

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