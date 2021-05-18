function [stat,x,y,z,Re,fltype,dstar,rlam,spanv,sumw]=read3Dxystats(filename)
%
% [stat,x,y,z,Re,fltype,dstar,rlam,spanv,sumw]=read3Dxystats(filename)
% 
% Read 3Dxy-statistics file as defined in Simson.
% See wxys.f and boxxys.f for details.
%
% This version uses the vector form of fread
%
% The output "stat" is a structured object with
% all the statistics. (NOTE: function under construction).
%
%NOTE: boxheight is along the y-axis.
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
bau=fread(fid,1,'int');
xl=fread(fid,1,'float64');
zl=fread(fid,1,'float64');
t=fread(fid,1,'float64');
shift=fread(fid,1,'float64');eol=fread(fid,2,'int');
A=fread(fid,1,'uint8=>char');
mhd_n=fread(fid,1,'float64');
b0=fread(fid,3,'float64');eol=fread(fid,2,'int');
nx=fread(fid,1,'int');
nyp=fread(fid,1,'int');
nzc=fread(fid,1,'int');
nfzsym=fread(fid,1,'int'); eol=fread(fid,2,'int');
fltype=fread(fid,1,'int');
dstar=fread(fid,1,'float64'); eol=fread(fid,2,'int');
if fltype<0
     rlam=fread(fid,1,'float64'); eol=fread(fid,2,'int');
elseif fltype >=6
    bstart=fread(fid,1,'float64');
    blenght=fread(fid,1,'float64');
    rlam=fread(fid,1,'float64');
    spanv=fread(fid,1,'float64'); eol=fread(fid,2,'int');
end
sumw=fread(fid,1,'float64');
nxys=fread(fid,1,'int'); eol=fread(fid,1,'int');
nz=nzc/(1-nfzsym);

% List of velocity statistics.
Lvs={'u','v','w','u2','v2','w2','omX','omY','omZ','omX2','omY2','omZ2', ...
    'uv','uw','vw'};
last=6; %nxys

for i=1:1:last
    eol=fread(fid,1,'int');
    dummy=fread(fid,nx*nyp*nz,'float64');
    vel = fliplr(reshape(dummy,nx,nyp,nz)); 
    ix = [nx/2:nx, 1:nx/2-1];
    vel(:,:,:) = vel(ix,:,:);
    stat.(Lvs{i})=vel;
    eol=fread(fid,1,'int');
end
fclose(fid);
Re=Re*dstar;
sumw=sumw/dstar;
x=(0:nx-1)/nx * xl/dstar;
y=(1-cos(linspace(0,pi,nyp))) * 1/2 * 2.0/dstar;
z= (-nz/2:1:nz/2-1)/nz * zl/dstar;

 end