function sigma = nmfm_psi1(coeffs, p0, p1, Deltas)
    sigma = 0;
    for i=1:size(coeffs,2)
        sigma = sigma + p1*Deltas{i}*coeffs(1:length(p0),i)/i;
    end
end
