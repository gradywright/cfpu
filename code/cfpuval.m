function [potential,X,Y,Z] = cfpuval(puminfo,gridsize)

% Local variables
x = puminfo.x;
y = puminfo.y;  M = size(y,1);
patchRad = puminfo.patchRad;
idxinfo = puminfo.idxinfo;
patchinfo = puminfo.patchinfo;

% Set up meshgrid for computing the potential
[minx,maxx] = bounds(x);

griddx = (maxx-minx)/gridsize;
griddx = max(griddx);
startx = (minx(1)-3*griddx);  endx = (maxx(1)+3*griddx);
starty = (minx(2)-3*griddx);  endy = (maxx(2)+3*griddx);
startz = (minx(3)-3*griddx);  endz = (maxx(3)+3*griddx);

xx = startx:griddx:endx;
yy = starty:griddx:endy;
zz = startz:griddx:endz;
[X, Y, Z] = meshgrid(xx,yy,zz);
% % Flatten vectors
% xe = [X(:) Y(:) Z(:)];
% m = length(xe);
[mmy,mmx,mmz] = size(X);
m = mmx*mmy*mmz;

% Determine which evaluation nodes belong to which patch
% [idxe,De] = rangesearch(xe,y,patchRad);
% [idxe,De] = determineGridsInPatches(y,X,Y,Z,patchRad,startx,dx);
% [eval_vec,temp_De] = rangesearch(tree,xe,patchRad);
% idxe = cell(M,1);
% % De = cell(M,1);
% % temp_length = cellfun(@length,eval_vec);
% % id_length = find(id > 0);
% for j=1:m
%     id = eval_vec{j};
%     for k = 1:length(id)
%         idxe{id(k)} = [idxe{id(k)} j];
% %         De{id(k)} = [De{id(k)} temp_De{j}];
%     end
% end

% Radial kernels to use
eta = puminfo.kernelinfo.eta;       % Curl-free 
phi = puminfo.kernelinfo.phi;       % Exact interpolation of potential
order = puminfo.kernelinfo.order;   % Curl-free polynomial degree

% Regularization:
exactinterp = puminfo.reginfo.exactinterp;

% Variables for sparse storage
% idxe_vec = [idxe{:}]';
% idxe_vec = vertcat(idxe{:});
idxe_patch = cell(1,M);
patch_vec = cell(1,M);

% Perform RBF locally for each patch and store
Psi = cell(1,M);              % values for pum weight function on each patch
potential_local = cell(1,M);  % values for potential function

% Take each column of the patches to aid in better parfor performance.
patchx = y(:,1);
patchy = y(:,2);
patchz = y(:,3);

% Loop over each patch and store local interpolant and pum weight function
parfor k = 1:M
%     warning('off','MATLAB:nearlySingularMatrix')
%     if mod(k,40) == 0
%         fprintf('k = %d out of M = %d\n',k,M)
%     end
    x_local = x(idxinfo(k).id,:); % Grab nodes on patch
    xx_local = x_local(:,1).';
    xy_local = x_local(:,2).';
    xz_local = x_local(:,3).';
    n = size(x_local,1);           % Number of nodes on patch
    
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
%     temp = De{k}.'/patchRad;
%     De = sqrt(sum((y(k,:)-xe_local).^2,2));
%     Psi{k} = (max(1-De/patchRad,0));
%     Psi{k} = (max(1-temp,0).^4).*(4*temp+1);
    Psi{k} = weight(De,patchRad,0);
%     Psi{k} = weight(De,patchRad,0);
    
    mlocalx = length(xx);
    mlocaly = length(yy);
    mlocalz = length(zz);
    xx = repmat(xx,[mlocaly 1 mlocalz]);
    yy = repmat(yy,[1 mlocalx mlocalz]);
    zz = repmat(zz,[mlocaly mlocalx 1]);
    xe_local = [xx(id) yy(id) zz(id)];
%     xe_local = xe(idxe{k},:);      % Evaluation points on each patch
    mm = size(xe_local,1);         % Number of evaluation points
    if mm == 0
        continue
    end
        
    % Make things faster by setting temp variables.
    coeffsp = patchinfo(k).coeffsp;
    coeffsx = patchinfo(k).coeffsx;
    coeffsy = patchinfo(k).coeffsy;
    coeffsz = patchinfo(k).coeffsz;
    coeffs_correction = patchinfo(k).coeffs_correction;
    coeffs_correction_const = coeffs_correction(end);
    coeffs_correction = coeffs_correction(1:end-1);

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
