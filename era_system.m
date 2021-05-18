function [A,B,C,D,Sig,Ad,Bd,Cd,Dd,U,V,H1] = era_system(yyin,dt,mc,mo,nr,varargin)
%
%   [A,B,C,D,Sig,Ad,Bd,Cd,Dd] = era_system(y,dt,mc,mo,nr[,toll,sampl])
%

%Input
if nargin > 5
    toll = varargin{1};
else
    toll = eps;
end

if nargin > 6
    sampl = varargin{2};
else
    sampl = 1;
end

nu = size(yyin,3); ny = size(yyin,1);

% downsampling
mc = floor(mc/sampl);
mo = floor(mo/sampl);
nt = floor((size(yyin,2)-1)/sampl)*sampl;%-1
yy = zeros(ny,length(1:sampl:nt)+1,nu); yy1 = yy;
for i = 1:sampl
    yy1(:,2:end,:) = yy1(:,2:end,:) + yyin(:,i+1:sampl:nt+1,:)/sampl;
    yy (:,2:end,:) = yy (:,2:end,:) + yyin(:,i  :sampl:nt+0,:)/sampl;
end

%Buildin H and H1
H = zeros(mo*ny,mc*nu); H1 = H;
[itc,ito] = meshgrid((1:mc),(1:mo));
for iy = 1:ny
    for iu = 1:nu
        y  = yy(iy,:,iu);
        y1 = yy1(iy,:,iu);
        I = iy+(0:mo-1)*ny;
        J = iu+(0:mc-1)*nu;
        H (I,J) = y ((ito+itc));
        H1(I,J) = y1((ito+itc));
    end
end

%SVD
[U,Sig,V] = svd(H,'econ'); Sig = diag(Sig);
%opts = [];
%[U,Sig,V] = lmsvd(H,rank(H),opts); Sig = diag(Sig);

%Reduction
nr = abs(min([nr,find(Sig/Sig(1)<=toll,1,'first')]))
Ur = U(:,1:nr); Vr = V(:,1:nr); sigr = Sig(1:nr);

SIG = spdiags(sqrt(sigr),0,nr,nr);
iSIG = spdiags(1./sqrt(sigr),0,nr,nr);

%Time discrete
Ad = iSIG*Ur'*H1*Vr*iSIG;
Bd = SIG*Vr'; Bd = Bd(:,1:nu);
Cd = Ur*SIG;  Cd = Cd(1:ny,:);
Dd = zeros(ny,nu);

%Time continue
A = real(logm(Ad))/dt; B = Bd; C = Cd; D = Dd;