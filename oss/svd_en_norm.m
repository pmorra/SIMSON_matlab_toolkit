function [U,S,V] = svd_en_norm(F,invF,H)

% Computes the svd w.r.t. the energy norm
%
% INPUT:  F, Cholesky factorization matrix for energy norm
%         invF, inverse of F
%         H, matrix to decompose with the svd
%
% OUTPUT: U, matrix of left singular vectors
%         S, matrix of singular values
%         V, matrix of right singular vectors
%
% Pierluigi Morra, 2020

% SVD of the resolvent w.r.t. the energy norm
[LSV,S,RSV] = svd(F*H*invF);

% Remove weights due to the energy norm
U = invF*LSV;
V = invF*RSV;

end