function [F,invF,W] = cheb_w_en_norm(ny,extr)

% Computes useful matrices related to the energy norm
% on a Gauss-Lobatto grid (along y)
%
% INPUT:  ny,   the number of points of the original U(y)
%         extr, if true includes the extrema at y(1) and y(ny)
%
% OUTPUT: W,    matrix of weights for the enrgy norm u'*W*u
%         F,    Cholesky factorization of W
%         invF, inverse of F
%
% Pierluigi Morra, 2020


[~,w] = clencurt(ny-1);

W = diag([w w w]);
F = chol(W);
invF = diag(1./diag(F));

if ~extr
  sel = [2:ny-1 ny+(2:ny-1) 2*ny+(2:ny-1)];
  W = W(sel,sel);
  F = F(sel,sel);
  invF = invF(sel,sel);
end

end