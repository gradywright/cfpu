idp = 2;
idg = round((y(idp,:)-(minx-3*dx))/dx)+1;
ix = idg(1); iy = idg(2); iz = idg(3);
%[X(ix,iy,iz) Y(ix,iy,iz) Z(ix,iy,iz)]
factor = round(patchRad/dx);
ixs = ix-factor:ix+factor;
iys = iy-factor:iy+factor;
izs = iz-factor:iz+factor;
plot3(y(idp,1),y(idp,2),y(idp,3),'b.','MarkerSize',12);
hold on
tempx = X(iys,ixs,izs);
tempy = Y(iys,ixs,izs);
tempz = Z(iys,ixs,izs);
xep = [tempx(:) tempy(:) tempz(:)];
plot3(xep(:,1),xep(:,2),xep(:,3),'k.')
plot3(X(iy,ix,iz),Y(iy,ix,iz),Z(iy,ix,iz),'cx')
xe_local = xe(idxe{idp},:);
plot3(xe_local(:,1),xe_local(:,2),xe_local(:,3),'m.','MarkerSize',10)
id = sum((y(idp,:)-xep).^2,2) < patchRad^2;
id2 = sub2ind(size(X),iys(id),ixs(id),izs(id));
xe_local2 = xe(id2,:);
plot3(xe_local2(:,1),xe_local2(:,2),xe_local2(:,3),'g.','MarkerSize',10)
hold off
