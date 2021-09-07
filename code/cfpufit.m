function puminfo = cfpufit(x,nrml,y,kernelinfo,reginfo)

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

% Matrices/vectors that are used over and over again
zm = zeros(l);
zv = zm(:,1);

% Array of structures for holding the patch data
puminfo = struct('x',x,'y',y,'patchRad',patchRad,'idxinfo',[],'patchinfo',[],'kernelinfo',kernelinfo,'reginfo',reginfo);
patchinfo = repmat(struct(...
    'n',0,'coeffsx',[],'coeffsy',[],'coeffsz',[],'coeffsCorrection',[]),M,1);
idxinfo = repmat(struct('id',[]),M,1);

opts.SYM = true;

% Loop over each patch and store local interpolant and Wendland function
parfor k = 1:M
%     if mod(k,40) == 0
%         fprintf('k = %d out of M = %d\n',k,M)
%     end
    id = idx{k};
    idxinfo(k).id = id;
    x_local = x(id,:);      % Grab nodes on patch
    xx_local = x_local(:,1).';
    xy_local = x_local(:,2).';
    xz_local = x_local(:,3).';

    n = size(x_local,1);        % Number of nodes on patch
    patchinfo(k).idx = idx{k};
    patchinfo(k).n = n;
    
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

    patchinfo(k).coeffsp = coeffsp;
    patchinfo(k).coeffsx = coeffsx;
    patchinfo(k).coeffsy = coeffsy;
    patchinfo(k).coeffsz = coeffsz;

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
    
    patchinfo(k).coeffs_correction = coeffs_correction;
end

puminfo.patchinfo = patchinfo;
puminfo.idxinfo = idxinfo;

end
