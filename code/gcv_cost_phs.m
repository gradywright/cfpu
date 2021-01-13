function score = gcv_cost_phs(lam,A,P,b,n)

Atilde = A + lam*eye(n);
temp = Atilde\P;
schurC = (P'*temp)\P';
S = (A - (A*temp - P)*schurC)/Atilde;
score = 1/n*sum((S*b-b).^2)./(1-trace(S)/n).^2;

end
