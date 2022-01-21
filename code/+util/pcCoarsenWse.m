function xc = pcCoarsenWse(x,M,area)
% pcCoarsenWse Coarsens a 3D point cloud using Weighted Sample Elimination (WSE)
% as implemented in Eliminate function of the cyCodeBase package.
%
% xc = pcCoarsenWse(x,M,area) coarsens the point cloud x (stored as a
% N-by-3 array of doubles) to exactly M points. area is an estimate for the 
% surface area of the object.
%
% xc = pcCoarsenPoissonDisk(x,M) same as above but area is approximated using
% the bounding box.

% Do some checking of the input parameters before calling the mex file

if nargin < 2
    error('The input point cloud X and the number of coarsened points M must be specified.');
end

if nargin == 2
    % Set the area to zero and let the cyCode estimate the value
    area = 0;
end

[N,d] = size(x);

if N <= 1 || d < 3
    error('The input point cloud should be an N-by-3 array, where N is the number of points, which should be greater than 1');
end

if M > N
    error('The number of coarsened points should be less than the number of points in the point cloud')
end

xc = util.mex_pcCoarsenWse(x,M,area);

end

