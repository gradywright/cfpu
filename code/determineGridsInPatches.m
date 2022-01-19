function [ide_patch,dist_patch] = determineGridsInPatches(patches,X,Y,Z,patchRad,startx,dx)

[mmy,mmx,mmz] = size(X);

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

    r2 = (patchx(j)-xx).^2 + (patchy(j)-yy).^2 + (patchz(j)-zz).^2;
    id = r2(:) < patchRad^2;
    
    ixs2 = repmat(ixs,[length(yy) 1 length(zz)]);
    iys2 = repmat(iys.',[1 length(xx) length(zz)]);
    izs2 = repmat(shiftdim(izs,-1),[length(yy) length(xx) 1]);
    
    temp_idg = (iys2 + (ixs2 - 1).*mmy) + (izs2-1)*(mmx*mmy);
    
    temp_idg = temp_idg(id);
    r2 = sqrt(r2(id));
    id = temp_idg > 0;
    temp_idg = temp_idg(id);
    r2 = r2(id);
    ide_patch{j} = temp_idg;
    dist_patch{j} = r2;
end

