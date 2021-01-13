function phi = weight(r,delta,k)

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

% b = 2;
% r2 = r.^2;
% 
% if k == 0
%    phi = exp(-b*r2./max(delta^2-r2,0));
% elseif k == 1
%    denom = max(delta^2-r2,0);
% %    phi = (-2*delta^2*b)*(r.*exp(-b*r2./denom)./denom.^2);
%    phi = (-2*delta^2*b)*(exp(-b*r2./denom)./denom.^2);
% else
%    error('PU Weight function error');
% end
   
end