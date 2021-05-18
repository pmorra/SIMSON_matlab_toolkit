% ***********************************************************************
%
% $HeadURL: https://www2.mech.kth.se:/svn/simson/trunk/matlab/readdns.m $
% $LastChangedDate: 2009-04-29 13:51:52 +0200 (Wed, 29 Apr 2009) $
% $LastChangedBy: pschlatt@MECH.KTH.SE $
% $LastChangedRevision: 1404 $
%
% ***********************************************************************
function [vel,x,y,z,Lx,Ly,Lz,fr,Re,flowtype,dstar]=readfou(filename)
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

%
% Define parameters based on retreived data
%
nx=storl(1);
ny=storl(2);
nz=storl(3);
Re=Re;
Lx=length(1);
Lz=length(2);
Ly=2/dstar;           % In bls, dstar is defined as 2/boxsize
fr=length(3);

realpos=kron(ones(1,nx),[1 0])==1;

vel.u = zeros(nx,ny,nz) + 1i*zeros(nx,ny,nz);
vel.v = zeros(nx,ny,nz) + 1i*zeros(nx,ny,nz);
vel.w = zeros(nx,ny,nz) + 1i*zeros(nx,ny,nz);

disp(' - Reading u');
for indz=1:nz
  for indy=ny:-1:1
    fread(fid,1,'int');
    vec=fread(fid,nx*2,'float64');
    vel.u(:,indy,indz) = vec(realpos) + 1i*vec(~realpos);
    fread(fid,1,'int');
  end
end

disp(' - Reading v');
for indz=1:nz
  for indy=ny:-1:1
    fread(fid,1,'int');
    vec=fread(fid,nx*2,'float64');
    vel.v(:,indy,indz) = vec(realpos) + 1i*vec(~realpos);
    fread(fid,1,'int');
  end
end

disp(' - Reading w');
for indz=1:nz
  for indy=ny:-1:1
    fread(fid,1,'int');
    vec=fread(fid,nx*2,'float64');
    vel.w(:,indy,indz) = vec(realpos) + 1i*vec(~realpos);
    fread(fid,1,'int');
  end
end
fclose(fid);

scale=1/dstar;
x=Lx/nx*(-nx/2:1:nx/2-1)';
if nz > 1
    z=Lz/nz*(-floor(nz/2):1:floor(nz/2)-1)';
else
    z=0;
end
y=scale*(1-cos(pi*(0:1/(ny-1):1)))';
x=-x(1)+x;