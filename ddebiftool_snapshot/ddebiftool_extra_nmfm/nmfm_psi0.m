function sigma = nmfm_psi0(coeffs, p0, p1, Deltas)
    sigma = 0;
    for i=1:size(coeffs,2)
        sigma = sigma + p0*Deltas{i}*coeffs(1:length(p0),i)/i + p1*Deltas{i+1}*coeffs(1:length(p0),i)/(i*(i+1));
    end
end
