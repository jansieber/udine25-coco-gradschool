function d=dde_psol_mult_dist(x,y)
%% log distance between two potential Floquewt multipliers in the complex plane
% distance between 0.1 and 1 should be greater than between 1 and 2, so for
% measuring the distance between Floquet multipliers we take the logarithm
% for the imaginary part we wrap around at pi instead of 0
xl=log(x);
yl=log(y);
imw=@(z)1i*(mod(imag(z)+pi,2*pi)-pi);
xl=real(xl)+imw(xl);
yl=real(yl)+imw(yl);
d=abs(xl-yl);
end
