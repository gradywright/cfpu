function [CP,P] = curlfreePoly(x,l)
% curlfreePoly Basis for curl-free polynomials up to degree l evaluated at x
%
% [CP,P] = curlfreePoly(x,l) Evaluates a basis for curl-free polynomials up to
% degree l at the values in x.  Returns the results in CP, which has 3*n rows
% and (l+1)*(l+2)*(l+2)/6-1 columns.  The kth 3-tuple of rows contains the
% curl-free polynomial evaluated at the kth value of x.  P contains the
% potentials for each curl-free polynomial evaluated at x.
%
% Note only degrees l = 1 and 2 are currently implemented.

% Copyright 2022 by Grady B. Wright

n = size(x,1);

if l == 1
    CP = zeros(3*n,3);
    P = zeros(n,3);
    P(:,1:3) = x;
    CP(:,1:3) = [repmat([1;0;0],n,1) repmat([0;1;0],n,1) repmat([0;0;1],n,1)];
elseif l == 2
    CP = zeros(3*n,9);
    P = zeros(n,9);
    P(:,1:3) = x;
    CP(:,1:3) = [repmat([1;0;0],n,1) repmat([0;1;0],n,1) repmat([0;0;1],n,1)];

    P(:,4:6) = 0.5*x.^2;
    zv = zeros(1,n);
    temp = [x(:,1)';zv;zv];
    CP(:,4) = temp(:);
    temp = [zv;x(:,2)';zv];
    CP(:,5) = temp(:);
    temp = [zv;zv;x(:,3)'];
    CP(:,6) = temp(:);

    P(:,7) = x(:,2).*x(:,3);
    temp = [zv;x(:,3)';x(:,2)'];
    CP(:,7) = temp(:);

    P(:,8) = x(:,1).*x(:,3);
    temp = [x(:,3)';zv;x(:,1)'];
    CP(:,8) = temp(:);

    P(:,9) = x(:,1).*x(:,2);
    temp = [x(:,2)';x(:,1)';zv];
    CP(:,9) = temp(:);
else
    error('Degree not implemented')
end

end