function [Rex,Re_delta,xa,delta,x_0] = calc_rex(x,Re)
%
%   [Rex,Re_delta,xa,delta,x_0] = calc_rex(x,Re)
%

d1_0  = 1;
U     = 1;
nu    = d1_0 * U / Re;

x_0   = - (d1_0/1.72)^2 * U / nu;
xa    = x - x_0;

Rex   = xa * U / nu;

delta = d1_0 * sqrt(xa / xa(1));
Re_delta = delta * U / nu;

end