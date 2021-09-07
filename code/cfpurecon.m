function [potential,X,Y,Z] = cfpurecon(x,nrml,y,kernelinfo,reginfo,gridsize)

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

% Fitting parameters:
regularization = reginfo.regularization;
lambda = reginfo.lambda;
schurcmplmnt = reginfo.schurcmplmnt;
exactinterp = reginfo.exactinterp;
regularizationi = reginfo.regularizationi;
lambdai = reginfo.lambdai;

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
zm = zeros(l);
zv = zm(:,1);

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
%     if mod(k,40) == 0
%         fprintf('k = %d out of M = %d\n',k,M)
%     end
    id = idx{k};
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
    if regularization ~= 3
        % Handle any potential singularities in zeta
        zeta_temp(1:n+1:end) = 0;
    else
        % Handle any potential singularities in zeta
        zeta_temp(r == 0) = 0;
    end
    
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
    
    if regularization == 1    % gcv regularization
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
        [lam,~,flag,~] = fminbnd(@gcv_cost_phs3,-10,35,[],z,D,3*n);
        lam = 3*n*exp(-lam);

        A = A + lam*eye(3*n);
        
        % Compute the coefficients
        coeffs = F2*(U*(z./(D.^2 + lam)));
        coeffsp = G1\(w1-F1'*(A*coeffs));       
    else
        if regularization == 2   % Manual regularization
            A(1:3*n,1:3*n) = A(1:3*n,1:3*n) + 3*n*lambda*eye(3*n);
        elseif regularization == 4
            if any(trbl_id(idx{k}))
                A(1:3*n,1:3*n) = A(1:3*n,1:3*n) + 3*n*lambda*eye(3*n);
            end
        end
        % Compute the coefficients
        if schurcmplmnt == 0
%             coeffs = [A CFP;CFPt zm]\b;
%             coeffs = A\b;
            coeffs = linsolve(A,b,opts);
            coeffsp = coeffs(3*n+1:end);
            coeffs = coeffs(1:3*n);        
        else
            A = A(1:3*n,1:3*n);
            b = b(1:3*n);
            coeffsp = pinv(CFPt*(A\CFP))*(CFPt*(A\b));
            coeffs = A\(b-CFP*coeffsp);
        end
    end
    
    % Make things faster by setting temp variables.
    coeffsx = coeffs(1:3:3*n).';
    coeffsy = coeffs(2:3:3*n).';
    coeffsz = coeffs(3:3:3*n).';    

    temp_potential_nodes = sum(eta_temp.*(dx.*coeffsx + dy.*coeffsy + dz.*coeffsz),2) + P*coeffsp;
    
    if exactinterp
        P = ones(n,1);
        A = ones(n+1);
        A(1:n,1:n) = phi(r);
        A(end) = 0;
        b = [temp_potential_nodes;0];
        if regularizationi == 1
            L = size(P,2);
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
            [lam,~,flag,~] = fminbnd(@gcv_cost_phs3,-10,35,[],z,D,n);
            lam = n*exp(-lam);
            A = A + lam*eye(n);
            % Compute the coefficients
            temp = F2*(U*(z./(D.^2 + lam)));
            coeffs_correction = [temp;G1\(w1-F1'*(A*temp))];
        else
            if regularizationi == 2
                A = A + n*lambdai*eye(n);
            elseif regularizationi == 4
                if any(trbl_id(idx{k}))
                    A = A + n*lambdai*eye(n);
                end
            end
%             coeffs_correction = [A P;P.' 0]\[temp_potential_nodes;0];
%             coeffs_correction = A\b;
            coeffs_correction = linsolve(A,b,opts);
        end
    else
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
