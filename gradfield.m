% ***********************************************************************
%
% $HeadURL: https://www2.mech.kth.se/svn/simson/trunk/matlab/gradfield.m $
% $LastChangedDate: 2006-11-16 21:05:30 +0100 (Thu, 16 Nov 2006) $
% $LastChangedBy: mattias@MECH.KTH.SE $
% $LastChangedRevision: 336 $
%
% ***********************************************************************
function [dxfou,dyfou,dzfou]=gradfield(fou,yF,kxvec,kzvec,varargin);
%
% Differentiate field in Fourier space
% Return data in Fourier space
%
% fou Two-dimensional plane in Fourier space
%

% Parse varargin
ddx = true; ddy = true; ddz = true;
if nargin == 5
  if length(varargin{1}) ~= 3
    ERRMSG = [' Error: the input to define the derivative to be computed',...
              ' must be an array with three entries.'];
    error(ERRMSG)
  else
    ddx = varargin{1}(1); ddy = varargin{1}(2); ddz = varargin{1}(3);
  end
elseif nargin > 5
    ERRMSG = 'Error: too many input arguments';
    error(ERRMSG);
end

s=size(fou);
Nx=2*s(1);
Nz=s(2)+1;

if length(s)>=3; Ny=s(3); else; Ny=1;end;
if length(s)>=4; Ns=s(4); else; Ns=1;end;
if length(s)>=5; Nr=s(5); else; Nr=1;end;

NNy=Ny/3;

%
% Set up differentiation matrix including boundaries
%
[xx,DM]=chebdif(NNy,1);

dxfou=0.0*fou;
dyfou=0.0*fou;
dzfou=0.0*fou;

%
% Compute streamwise and spanwise derivatives
%
for indy=1:NNy
  if ddz
    for indz=1:Nz-1
    % d/dz
    dzfou(:,indz,indy)       = sqrt(-1)*kzvec(indz)*fou(:,indz,indy);
    dzfou(:,indz,indy+NNy)   = sqrt(-1)*kzvec(indz)*fou(:,indz,indy+NNy);
    dzfou(:,indz,indy+2*NNy) = sqrt(-1)*kzvec(indz)*fou(:,indz,indy+2*NNy);
    end
  end
  if ddx
    for indx=1:Nx/2
    % d/dx
    dxfou(indx,:,indy)       = sqrt(-1)*kxvec(indx)*fou(indx,:,indy);
    dxfou(indx,:,indy+NNy)   = sqrt(-1)*kxvec(indx)*fou(indx,:,indy+NNy);
    dxfou(indx,:,indy+2*NNy) = sqrt(-1)*kxvec(indx)*fou(indx,:,indy+2*NNy);
    end
  end
end

%
% Wall-normal derivatives are taken in physical space
%
if ddy
  % 3/2 rule
  Nx_extra = 0;%Nx*0.5;
  Nz_extra = 0;%Nz*0.5;

  foup=fou2phys(fou,Nx_extra,Nz_extra);
  for indx=1:Nx+Nx_extra
    for indz=1:Nz+Nz_extra
      tmpfou(1:NNy)=foup(indx,indz,1:NNy);
      dyfoup(indx,indz,1:NNy)=(DM*tmpfou');

      tmpfou(1:NNy)=foup(indx,indz,NNy+1:2*NNy);
      dyfoup(indx,indz,NNy+1:2*NNy)=(DM*tmpfou');

      tmpfou(1:NNy)=foup(indx,indz,2*NNy+1:3*NNy);
      dyfoup(indx,indz,2*NNy+1:3*NNy)=(DM*tmpfou');
    end
  end
  %dyfou=phys2fou(dyfoup);
  dummy = phys2fou(dyfoup); %size(dummy)
  dyfou = dummy(1:Nx/2,[1:Nz/2 Nz/2+Nz_extra+1:Nz+Nz_extra-1],:);
end
%
% Compute wall-normal derivatives in Fourier space
% Does not work at the moment???
%
%for indx=1:Nx/2
%  for indz=1:Nz-1
    % d/dy

%    tmpfou(1:NNy)=fou(indx,indz,1:NNy);
%    dyfou(indx,indz,1:NNy)=(DM*tmpfou');
    
%    tmpfou(1:NNy)=fou(indx,indz,NNy+1:2*NNy);
%    dyfou(indx,indz,NNy+1:2*NNy)=(DM*tmpfou');

%    tmpfou(1:NNy)=fou(indx,indz,2*NNy+1:3*NNy);
%    dyfou(indx,indz,2*NNy+1:3*NNy)=(DM*tmpfou');
%  end
%end
