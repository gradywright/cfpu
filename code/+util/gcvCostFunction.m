function score = gcvCostFunction(lam,z,d,n)
% gcvCostFunction Cost function for determining the optimal smoothing spline
% parameter based on generalized cross validation.

% Copyright 2022 by Grady B. Wright

lam = exp(-lam);
temp = (n*lam)./(d.^2 + n*lam);
score = n*sum((temp.*z).^2)./sum(temp).^2;

end
