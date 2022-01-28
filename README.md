# CFPU
Curl-Free RBF Partition of Unity (**_CFPU_**) method for implicit surface reconstruction from oriented point clouds
![ScreenShot](https://raw.github.com/gradywright/cfpu/master/ptcloud_ex.png)![ScreenShot](https://raw.github.com/gradywright/cfpu/master/cfpu_ex.png)

This is a MATLAB implementation of the CFPU method described in [1].  The code can be used to produce an implicit (zero level-set) surface to an oriented point cloud - a point cloud consisting of points and (approximate) normals.

The following example shows how to use the code to produce the surface displayed above:
```
% Load the point cloud
load('../ptclouds/homer.mat');
% Reconstruct surface using CFPU
[potential,X,Y,Z] = cfpurecon(ptcloud,normals,patches,256);
% Plot the surface
fv = isosurface(X,Y,Z,potential,0); p = patch(fv); isonormals(X,Y,Z,potential,p);
set(p,'FaceColor',0.9*[1 1 1],'EdgeColor','none')
daspect([1 1 1]), lighting phong, material dull, view([90 0]);
camlight('right','infinite'), axis off tight, set(gcf,'Color',[1 1 1]) 
```

For more examples from [1] see [cpfu_ex.m](https://raw.github.com/gradywright/cfpu/master/cpfu_ex.m)

Note that the examples use point cloud sets that were obtained from other sites, but stored as .mat files for convenience.  The Stanford Bunny, Happy Buddha, Stanford Dragon, and Armadillo data were obtained from the [Stanford University 3D Scanning Repository](http://graphics.stanford.edu/data/3Dscanrep/). The Homer, Raptor, Filigree, Pump Carter, Dancing Children, Gargoyle, and De Bozbezbozzel data were obtained from the [AIM@SHAPE-VISIONAIR Shape Repository](http://visionair.ge.imati.cnr.it).

The code also include MEX interface files to some useful point cloud utility functions from the [VCGLib Package](https://github.com/cnr-isti-vclab/vcglib) and [cyCodeBase Package](https://github.com/cemyuksel/cyCodeBase) that can be used in conjuction with the CFPU method.  You will need to download these libraries and compile the MEX files according to the instructions in the help text of the C++ mex files.

References:

[1] K. P. Drake, E. J. Fuselier, and G. B. Wright. Implicit Surface Reconstruction with a Curl-free Radial Basis Function Partition of Unity Method. [arXiv:2101.05940](https://arxiv.org/abs/2101.05940)
