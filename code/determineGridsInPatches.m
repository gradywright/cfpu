function [ide_patch,dist_patch] = determineGridsInPatches(patches,X,Y,Z,patchRad,startx,dx)

[mmy,mmx,mmz] = size(X);
% id_entire_grid = reshape(1:mmy*mmx*mmz,[mmy mmx mmz]);

% id_entire_grid2 = 0*id_entire_grid;
% for k = 1:length(x)
%     idx = round((x(k,:)-startx)/dx)+1;
%     ix = idx(1); iy = idx(2); iz = idx(3);
%     factor = 2;
%     ixs = max(ix-factor,1):min(ix+factor,mmx);
%     iys = max(iy-factor,1):min(iy+factor,mmy);
%     izs = max(iz-factor,1):min(iz+factor,mmz);
%     id_entire_grid2(iys,ixs,izs) = id_entire_grid(iys,ixs,izs);
% end
% id_entire_grid = id_entire_grid2;

ide_patch = cell(length(patches),1);
dist_patch = ide_patch;
patchx = patches(:,1);
patchy = patches(:,2);
patchz = patches(:,3);
parfor j = 1:length(patches)
    idg = round((patches(j,:)-startx)/dx)+1;
    ix = idg(1); iy = idg(2); iz = idg(3);
    factor = round(patchRad/dx);
    ixs = max(ix-factor,1):min(ix+factor,mmx);
    iys = max(iy-factor,1):min(iy+factor,mmy);
    izs = max(iz-factor,1):min(iz+factor,mmz);
    xx = (startx(1)+(ixs-1)*dx);
    yy = (startx(2)+(iys-1)*dx).';
    zz = shiftdim(startx(3)+(izs-1)*dx,-1);
%     err = [xx-X(1,ixs,1) yy-Y(iys,1,1)' zz-squeeze(Z(1,1,izs))'];
%     if norm(err,inf) > 100*eps
%         'foo'
%     end
%     r22 = (patches(j,1)-X(iys,ixs,izs)).^2 + (patches(j,2)-Y(iys,ixs,izs)).^2 + (patches(j,3)-Z(iys,ixs,izs)).^2;
    r2 = (patchx(j)-xx).^2 + (patchy(j)-yy).^2 + (patchz(j)-zz).^2;
    
%     if norm(r2(:)-r22(:),inf) > 100*eps
%         'foo'
%     end
    id = r2(:) < patchRad^2;
%     tempx = X(iys,ixs,izs);
%     tempy = Y(iys,ixs,izs);
%     tempz = Z(iys,ixs,izs);
%     xep = [tempx(:) tempy(:) tempz(:)];
    
    ixs2 = repmat(ixs,[length(yy) 1 length(zz)]);
    iys2 = repmat(iys.',[1 length(xx) length(zz)]);
    izs2 = repmat(shiftdim(izs,-1),[length(yy) length(xx) 1]);
    
%     temp_idg3 = sub2ind([mmy mmx mmz],ixs2,iys2,izs2);
    temp_idg = (iys2 + (ixs2 - 1).*mmy) + (izs2-1)*(mmx*mmy);
%     temp_idg = id_entire_grid(iys,ixs,izs);
    
    temp_idg = temp_idg(id);
    r2 = sqrt(r2(id));
    id = temp_idg > 0;
    temp_idg = temp_idg(id);
    r2 = r2(id);
%     id = sum((y(idp,:)-xep).^2,2) < patchRad^2;
    ide_patch{j} = temp_idg;
    dist_patch{j} = r2;
%     xep = [tempx(:) tempy(:) tempz(:)];
%     xx_reduced_grid = [xx_reduced_grid;xep];
end

