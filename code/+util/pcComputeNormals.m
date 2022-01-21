function nrmls = pcComputeNormals(x,m,k)
% pcComputeNormals Computes approximate normals for a 3D point cloud of a
% surface as implemented in the PointCloudNormal function of the VCGLIB package.
%
% nrmls = pcComputeNormals(x) computes approximate normals for
% the point cloud x (stored as a N-by-3 array of doubles) using 10 neighbors for
% each point and no smoothing sweeps.
%
% nrmls = pcComputeNormals(x,m,k) computes approximate normals for
% the point cloud x (stored as a N-by-3 array of doubles) using m neighbors for
% each point and k smoothing sweeps.

% Do some checking of the input parameters before calling the mex file
if nargin < 2
    m = 10;
    k = 0;
end

if nargin == 2
    k = 0;
end

[N,d] = size(x);

if N <= 1 || d < 3
    error('The input point cloud should be an N-by-3 array, where N is the number of points, which should be greater than 1');
end

if k < 0
    warning('Smoothing steps should be >= 0.  Setting this value to zero.');
    k = 0;
end

if m <= 2
    warning('Number of neighbors to use in computation should be > 2.  Setting this value to 10.');
    m = 10;
end
    
nrmls = util.mex_pcComputeNormals(x,m,k);

end

