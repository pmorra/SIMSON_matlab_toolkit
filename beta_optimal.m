function betaopt = beta_optimal(Re,x)

% Computes the beta that gives the max energy growth at x
% according to SIMSON non-dimensionalization
%
% Reference: Andersson et al. (1999), Physics of Fluids (11),134
  xa = x + Re/1.72^2;
  betaopt = sqrt(Re./xa)*0.45;
end