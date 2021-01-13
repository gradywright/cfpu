load HomerCFPUMPHSOrder2
fv = isosurface(X,Y,Z,g,0);
p = patch(fv);
isonormals(X,Y,Z,g,p);
daspect([1 1 1])
view(3)
% set(p,'FaceColor',[255 241 0]/255,'EdgeColor','none')
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
% set(p,'SpecularColorReflectance', 0, 'SpecularExponent', 100, 'DiffuseStrength', 0.8);
lighting phong
material dull
% view([180 -90])
% view([0 90])
% view([60 15])    % Dragon
% view([-90 4]);    % Armadillo
% view([90 0]);   % Buddha 
% view([-128 36]);    % cantius_tooth
view([90 0]);     % Children & Knot & Filigree
% view([0 90])
% view([93 4]);     % Bunny
% view([67 5]);     % Carter
% camlight headlight
camlight('right','infinite')
% camlight('headlight','infinite');  % Bunny
% l = light('Position',[100 10 0],'Style','infinite');
% l2 = light('Position',[-100 0 0],'Style','infinite');
axis off
set(gcf,'Color',[1 1 1]*1)
% view([-40 20])
% view([-30 10]);
axis tight
ax = axis;


%%
plot3(xc(:,1),xc(:,2),xc(:,3),'k.','MarkerSize',0.1)
view([0 90])
daspect([1 1 1])
axis off
shg
