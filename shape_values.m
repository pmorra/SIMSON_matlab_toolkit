function [ dstar,theta,H ] = shape_values(y,vel)
%shapeValues computes displacement and momentum thickness and shape factor
%   values along x direction.
%
% [ dstar,theta,H ] = shapeValues(y,vel)
%   y: as it comes from readdns
%  vel: it can be vel.u, where vel comes from mat2vel

dstar=trapz(y,(1-vel),2);
theta=trapz(y,(vel.*(1-vel)),2);
H=dstar./theta;
end
