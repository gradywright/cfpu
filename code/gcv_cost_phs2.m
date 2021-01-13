function score = gcv_cost_phs2(lam,ZQ,AZQ,C,D,b,n)

% Atilde = A + lam*eye(n);
% temp = Atilde\P;
% schurC = (P'*temp)\P';
% S = (A - (A*temp - P)*schurC)/Atilde;

temp1 = (ZQ'./(D+exp(-lam)));
temp2 = AZQ*temp1;
temp1 = ZQ*temp1;
S = temp2 + C*(eye(n) - temp2 - exp(-lam)*temp1);
score = 1/n*sum((S*b-b).^2)./(1-trace(S)/n).^2;

end
