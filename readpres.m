% ***********************************************************************
%
% $HeadURL: https://www2.mech.kth.se:/svn/simson/trunk/matlab/readdns.m $
% $LastChangedDate: 2009-04-29 13:51:52 +0200 (Wed, 29 Apr 2009) $
% $LastChangedBy: pschlatt@MECH.KTH.SE $
% $LastChangedRevision: 1404 $
%
% ***********************************************************************
function [pres,xF,yF,zF,Lx,Ly,Lz,t,Re,flowtype,dstar,pou,rlam,spanv]=readpres(filename)
%
% Read a velocity field in Fourier space as defined in Simson
% in variables u,v,w.
%
% This version uses the vector form of fread
%
% The velocity components u, v, w are concatenated
% on the third dimension of the output vel
%
% The odd ball is removed
%
rlam=0.0;
spanv=0.0;

%
% Open file
%
fid=fopen(filename,'r','ieee-le.l64');
eol=fread(fid,1,'int');

%
% If record is not 44 the file is either corrput or on big endian format
%
if eol ~= 44
  fclose(fid);
  disp(' ')
  disp(['Reading ' filename ' on big endian format'])
  fid=fopen(filename,'r','ieee-be.l64');
  eol=fread(fid,1,'int');
else
  disp(' ')
  disp(['Reading ' filename ' on little endian format'])
end
Re=fread(fid,1,'float64');
pou=fread(fid,1,'int');
length=fread(fid,4,'float64');
eol=fread(fid,2,'int');
storl=fread(fid,4,'int');
eol=fread(fid,2,'int');
flowtype=fread(fid,1,'int');
dstar=fread(fid,1,'float64');
eol=fread(fid,1,'int');
if flowtype==-1
  eol=fread(fid,1,'int');
  rlam=fread(fid,1,'float64');
  eol=fread(fid,1,'int');
elseif flowtype==-2
  eol=fread(fid,1,'int');
  rlam=fread(fid,1,'float64');
  spanv=fread(fid,1,'float64');
  eol=fread(fid,1,'int');
%elseif flowtype==4 || flowtype==5 || flowtype==6
%  eol=fread(fid,1,'int')
%  boundparam=fread(fid,2,'float64');
%  eol=fread(fid,1,'int')
%  bstart=boundparam(1)
%  blength=boundparam(2)
elseif flowtype>=4
  eol=fread(fid,1,'int');
  boundparam=fread(fid,4,'float64');
  eol=fread(fid,1,'int');
  bstart=boundparam(1);
  blength=boundparam(2);
  rlam=boundparam(3);
  spanv=boundparam(4);
end

%
% Define parameters based on retreived data
%
NNx=storl(1)/2;
NNy=storl(2);
NNz=storl(3);
Re=Re*dstar;
Lx=length(1)/dstar;
Lz=length(2)/dstar;
Ly=2/dstar;           % In bls, dstar is defined as 2/boxsize
t=length(3)/dstar;

realpos=kron(ones(1,NNx),[1 0]);

disp(' - Reading p');
for indz=1:NNz
  for indy=1:NNy
    fread(fid,1,'int');
    vec=fread(fid,NNx*2,'float64');
    rlu(:,indz,indy)=vec(~~realpos);
    ilu(:,indz,indy)=vec(~realpos);
    fread(fid,1,'int');
  end
end
fclose(fid);

scale=1/dstar;
padx=0;
padz=0;
NxF=2*NNx;
NzF=NNz;
xF=Lx/NxF*(-NxF/2:1:NxF/2)';
if NzF > 1
    zF=Lz/NzF*(-floor(NzF/2):1:floor(NzF/2))';
else
    zF=0;
end
yF=scale*(1+cos(pi*(0:1/(NNy-1):1)))';
xF=-xF(1)+xF;

%
% Shift velocity field in the streamwise direction in order
% to move the fringe to the end of the domain for spatial
% flows
%
kxvec=linspace(0,2*pi/Lx*(NNx-1),NNx);
if NNz > 1
    kzvec=linspace(0,2*pi/Lz*(NNz/2-1),NNz/2);
    kzvec=[kzvec -2*pi*NNz/2/Lz -fliplr(kzvec(2:end))];
else
    kzvec = 0;
end

xs = Lx/2.;
zs = 0.;

for i=1:NNx
  argx = -xs*kxvec(i);
  cx(i) = cos(argx);
  sx(i) = sin(argx);
end
for k=1:NNz
  argz = -zs*kzvec(k);
  for i=1:NNx
    ca(i)=cx(i)*cos(argz)-sx(i)*sin(argz);
    sa(i)=cx(i)*sin(argz)+sx(i)*cos(argz);
  end
  for j=1:NNy
    for i=1:NNx
      hr=rlu(i,k,j)*ca(i)-ilu(i,k,j)*sa(i);
      ilu(i,k,j)=ilu(i,k,j)*ca(i)+rlu(i,k,j)*sa(i);
      rlu(i,k,j)=hr;
    end
  end
end
%
%u=fftshift(u,1);
%v=fftshift(v,1);
%w=fftshift(w,1);

u=reshape(complex(rlu,ilu),NNx,NNz,NNy);


% Set to zero the zerozero mode for u component
%u(1,1,:)=u(1,1,:)*0.0;

%
% Concatenate the components on the 3rd dimension
%
pres=ccat(1,u);

%
% Remove odd ball
%
if NNz > 1
    pres=ccat(2,vel(:,1:NNz/2,:),vel(:,NNz/2+2:end,:));
%else
%    vel=ccat(2,vel(:,1,:));
end
