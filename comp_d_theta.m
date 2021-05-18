function [dstar,theta] = comp_d_theta(y,u)

% Assumes u(x,y);
ymax = max(abs(y));
ny = length(y);

[~,w] = clencurt(ny-1);
w = w*ymax/2;
dstar = ymax - u*w';
theta = u*w' - u.^2*w';
end