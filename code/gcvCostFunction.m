function score = gcv_cost_phs3(lam,z,d,n)

lam = exp(-lam);
temp = (n*lam)./(d.^2 + n*lam);
score = n*sum((temp.*z).^2)./sum(temp).^2;

end
