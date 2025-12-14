function h = nmfm_bt_create_solution(kappa,coeffs,Deltas,Dinv,devlt)

xi = kappa;
for i=1:size(coeffs,2)
    if size(coeffs,2) > 1
        coeffs(:,i) = coeffs(:,i)/i;
        xi = xi - Deltas{i}*coeffs(:,i);
     else % just a single vector
        xi = xi - Deltas{i}*coeffs;
    end
end
h = devlt([Dinv(xi),coeffs]);

end
