patchRad = puminfo.patchRad;
x = puminfo.x;
y = puminfo.y;
M = length(y);
plot3(x(:,1),x(:,2),x(:,3),'k.')
hold on
% plot3(y(:,1),y(:,2),y(:,3),'r.','MarkerSize',14)
for j=1:M
    [xx,yy,zz] = sphere(20);
    xx = xx*patchRad;
    yy = yy*patchRad;
    zz = zz*patchRad;
    xx = xx+y(j,1);
    yy = yy+y(j,2);
    zz = zz+y(j,3);
    h = surf(xx,yy,zz);
    set(h,'FaceAlpha',0.0,'FaceColor',0.8*[0 0 1],'EdgeColor',[66, 206, 245]/255)
end
plot3(x(:,1),x(:,2),x(:,3),'k.')
daspect([1 1 1])
view([90 0])
axis off

%%
pname = 'homer_patches_%03d';
for k = 1:360
    view([90+k 0])
    drawnow
end

%%
export_fig('-png','-m2',sprintf(pname,0));
set(gcf,'Color',[1 1 1]*1)
%%
for i = 1:180
   camorbit(2,0,'data',[0 0 1])
   set(gcf,'Color',[1 1 1]*1)
   drawnow
   export_fig('-png','-m2',sprintf(pname,i));
end