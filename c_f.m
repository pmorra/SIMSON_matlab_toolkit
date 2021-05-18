function [ c_f ] = c_f( y, arg, Re )
% c_f : computes the friction coefficient via Chebyshev polynomials
% cf = c_f(y,arg,Re)

addpath(pathdef);

ny = length(y);
Ly = y(end);

[~, D] = chebdif(ny,1);
Dw = D(1,:) * -2/Ly;

c_f = 2/Re*Dw*(arg);

end

