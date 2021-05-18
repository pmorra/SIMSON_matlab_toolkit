function [f] = frequency(N,spacelength)
% frequency: creates a 1D array for fft. The frequency are divided by 2*pi.
% INPUTS:
%           N :          total number of discrete points for the domain
%           spacelength: total lenght of the 1D domain which is being
%                        transformed into the fourier domain.
% OUTPUT:
%           f:           array of frequencies/(2*pi)
%
% Pierluigi Morra (2016-Nov)


f  = (-floor(N/2):-1+mod(N,2)+floor(N/2))/(spacelength);

end