function xc = pcCoarsenPoissonDisk(x,M)
% pcCoarsenPoissonDisk Coarsens a 3D point cloud using Poisson Disk Sampling as 
% implemented in PoissonPruning function of the VCGLIB package.
%
% xc = pcCoarsenPoissonDisk(x,M) coarsens the point cloud x (stored as a
% N-by-3 array of doubles) to approximately M points. Note that M is a very
% rough a approximation for the size of the coarsened set.  Typically, the
% coarsened point cloud has many more than M points.

% Do some checking of the input parameters before calling the mex file

if nargin < 2
    error('The input point cloud X and the number of coarsened points M must be specified.');
end

[N,d] = size(x);

if N <= 1 || d < 3
    error('The input point cloud should be an N-by-3 array, where N is the number of points, which should be greater than 1');
end

if M > N
    error('The number of coarsened points should be less than the number of points in the point cloud')
end

xc = util.mex_pcCoarsenPoissonDisk(x,M);

end

