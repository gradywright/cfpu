function phi = weight(r,delta,k)
% WEIGHT is the quadratic b-spline generator for the partition of unity  (PU)
% weights. This is combined with Shepard's method to construct the PU weights.
%
% phi = weight(r,delta,k) is the univariate quadratic b-spline function
% evaluated at the distance r, with support delta.  k=0 is the weight function
% and k=1 is the derivative of the weight function with respect to r.

% Copyright 2022 by Grady B. Wright

% Quadratic b-spline
r = r/delta;
phi = 0*r;

if k == 0
    id = r <= (1/3);
    phi(id) = 3/4 - 9/4*r(id).^2;
    id = (r > 1/3) & (r <= 1);
    phi(id) = 9/8*(1-r(id)).^2;
elseif k == 1
    id = r <= (1/3);
    phi(id) = -9/2/delta^2;
    id = (r > 1/3) & (r <= 1);
    phi(id) = (-9/4*(1-r(id))/delta^2).*(1./r(id));
else
   error('PU Weight function error');
end
   
end