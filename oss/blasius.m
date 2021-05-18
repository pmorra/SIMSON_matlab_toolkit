function [Y,U,Up,Upp] = blasius(Re,n,ymax)
%
% Computes the Blasius solution on Gauss-Lobatto points
%
% INPUT: Re,   Reynolds number (U_inf*delta0^*/nu)
%        n,    number of GL points
%        ymax, domain length
%
% OUTPUT: Y,   Gauss-Lobatto points
%         U,   velocity profile
%         Up,  first derivative dUdy
%         Upp, second derivative d2U/dy2
%
% Pierluigi Morra, 2020
%

% parameters
c0 = 1.7207876575715/sqrt(2);
c1 = 0.33205733621545*sqrt(2);

% a step size which is stable for runge-kutta
h0 = 0.05;

% Gauss-Lobatto points
for j = 1:n
  eta(j)=ymax*(1.-cos((j-1)*pi/(n-1)))/2;
end

% include parameters 
e = c0*eta;

% add points between GL points if step size > h0
g(1) = e(1);
nadd(1) = 1;
for k = 2:n
  de = e(k)-e(k-1);
  na = floor(de/h0)+1;
  ng = nadd(k-1);
  for l = 1:na
    g(ng+l) = e(k-1)+de*l/na;
  end
  nadd(k) = ng+na;
end

% B.C. at the wall for runge-kutta marching
x(1) = 0;
y(1) = 0;
z(1) = c1;
w(1) = 0;

% integration of the Blasius differential eq.
for i = 2:nadd(n)
  h = g(i)-g(i-1);
  [x(i),y(i),z(i)] = myrunge(x(i-1),y(i-1),z(i-1),h);
  w(i) = -x(i)*z(i)*c0*c0;
end

% allocate solutions U, Up, Upp
for k = 1:n
  u(k) = y(nadd(k));
  u1(k) = z(nadd(k))*c0;
  u2(k) = w(nadd(k));
end
Y = eta;
U = u;
Up = u1;
Upp = u2;

end


function [x1,y1,z1] = myrunge(x,y,z,h)
%
% Solves 1 step of the Blasius differential eq.
% with runge-kutta 4th order
%

f1 = @(x,y,z) (y);
f2 = @(x,y,z) (z);
f3 = @(x,y,z) (-x*z);

xt1 = h*f1(x,y,z);
yt1 = h*f2(x,y,z);
zt1 = h*f3(x,y,z);

xt2 = h*f1(x+xt1/2,y+yt1/2,z+zt1/2);
yt2 = h*f2(x+xt1/2,y+yt1/2,z+zt1/2);
zt2 = h*f3(x+xt1/2,y+yt1/2,z+zt1/2);

xt3 = h*f1(x+xt2/2,y+yt2/2,z+zt2/2);
yt3 = h*f2(x+xt2/2,y+yt2/2,z+zt2/2);
zt3 = h*f3(x+xt2/2,y+yt2/2,z+zt2/2);

xt4 = h*f1(x+xt3,y+yt3,z+zt3);
yt4 = h*f2(x+xt3,y+yt3,z+zt3);
zt4 = h*f3(x+xt3,y+yt3,z+zt3);

x1 = x + (xt1 +2*xt2+2*xt3+xt4)/6;
y1 = y + (yt1 +2*yt2+2*yt3+yt4)/6;
z1 = z + (zt1 +2*zt2+2*zt3+zt4)/6;

end
