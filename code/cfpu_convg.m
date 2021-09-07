% Smoothing parameter:
lambda = 1e-2;
regularization = 0;
exactinterp = 1;

%Define eta and zeta
phi = @(r) r;
eta = @(r) r.^3;
zeta = @(r) 3*r;  order = 2;
% eta = @(r) -r;
% zeta = @(r) -1./r; order = 1;


NN = [6144 8664 11616 18816 23064 27744 32856]';
load(sprintf('../ptclouds/trefoil_N%04d.mat',131424));
xe = xc;
m = length(xe);
load(sprintf('../ptclouds/trefoil_patches_M%d',864));
[~,nn_dist] = knnsearch(y,y,'k',2);
H = max(nn_dist(:,2));
delta = 3/4;
patchRad = (1 + delta)*H/2;
err = zeros(length(NN),2);
for nn = 1:length(NN)
    load(sprintf('../ptclouds/trefoil_N%04d.mat',NN(nn)));
    x = xc;
    N = size(x,1);
    nrml = nrml./(sqrt(sum(nrml.^2,2)));

    M = length(y);
    % Determine which nodes belong to which patch
    idx = rangesearch(x,y,patchRad);

    % Determine which evaluation nodes belong to which patch
    [idxe,De] = rangesearch(xe,y,patchRad);

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
    ns = zeros(M,1);
    % Loop over each patch and store local interpolant and Wendland function
    parfor k = 1:M
        x_local = x(idx{k},:); % Grab nodes on patch
        n = length(idx{k});    % Number of nodes on patch
        ns(k) = n;

        [CFP,P] = curlfreePoly(x_local,order);

        % Least squares
        if regularization == 3
            nlsq = min(max(ceil(n/5),20),n);
            nlsq = min(20,n);
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
        zeta_temp(r == 0) = 0;
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
            coeffs = [A CFP;CFPt zeros(l)]\[b;zeros(l,1)];
            coeffsp = coeffs(3*nlsq+1:end);
            coeffs = coeffs(1:3*nlsq);        
        end

        temp_potential_nodes = sum(eta_temp.*(dx.*coeffs(1:3:3*nlsq)' + dy.*coeffs(2:3:3*nlsq)' + dz.*coeffs(3:3:3*nlsq)'),2) + P*coeffsp;

        if exactinterp
            P = ones(n,1);
            dx = (x_local(:,1) - x_local(:,1).');
            dy = (x_local(:,2) - x_local(:,2).');
            dz = (x_local(:,3) - x_local(:,3).');
            r = sqrt(dx.^2 + dy.^2 + dz.^2);
            coeffs_correction = [phi(r) P;P.' zeros(1)]\[temp_potential_nodes;zeros(1,1)];
        else
%             P = [ones(n,1) P(:,1:3)];
            P = ones(n,1);
            coeffs_correction = P\temp_potential_nodes;
%             coeffs_correction = mean(temp_potential_nodes);
        end
        
        % Generate the distance matrix between the patch nodes and global evaluation points
        mm = length(idxe{k});
        xe_local = xe(idxe{k},:);

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
%                 potential_correction(idb) = [ones(mmb,1) P(:,1:3)]*coeffs_correction;
                potential_correction(idb) = ones(mmb,1)*coeffs_correction;
%                 potential_correction(idb) = potential_correction(idb) + coeffs_correction;
            end
        end

        % Calculate the weight function on each patch center and store it
        Psi{k} = weight(De{k}.',patchRad,0);
%         Psi{k} = (max(1-De{k}.'/patchRad,0));
        % Correct the potential function
        potential_local{k} = temp_potential - potential_correction;
%         potential_local{k} = temp_potential - mean(temp_potential_nodes);

        patch_vec{k} = k*ones(mm,1);
    end

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
    err(nn,1) = sqrt(mean(g.^2));
    err(nn,2) = norm(g,inf);
    fprintf('N=%05d, rms=%1.3e, max=%1.3e\n',N,err(nn,1),err(nn,2))
end


