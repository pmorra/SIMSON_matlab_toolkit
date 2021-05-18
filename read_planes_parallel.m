function [out,param] = read_planes_parallel(filename,maxNT);

% read_planes: reads the planes saved from SIMSON in parallel
%
% INPUT  : filename: name of the file
%          maxNT: maximum number of file to read
%
% OUTPUT : out:   cell array, in each cell is a plane at a time
%          param: axes (argx,argy), time series (ta), Reynolds (re) 
%
% Pierluigi Morra 2020
% email: pmorra@mech.kth.se
% 
tic
file_format = 'ieee-le.l64';

% Open the file (assuming it is little endian) 
fid = fopen(filename,'r',file_format);

% First line
eol = fread(fid,1,'int');
re = fread(fid,1,'float64');
pou = fread(fid,1,'int=>logical');
xl = fread(fid,1,'float64');
zl = fread(fid,1,'float64');
t = fread(fid,1,'float64');
xs = fread(fid,1,'float64');
eol = fread(fid,1,'int');

% Second line
eol = fread(fid,1,'int');
nxin = fread(fid,1,'int');
nypin = fread(fid,1,'int');
nzcin = fread(fid,1,'int');
nfzsin = fread(fid,1,'int');
eol = fread(fid,1,'int');

% Third line
eol = fread(fid,1,'int');
tpl = fread(fid,1,'int');
ivar = fread(fid,1,'int');
cpl = fread(fid,1,'float64');
fltype = fread(fid,1,'int');
dstar= fread(fid,1,'float64');
eol = fread(fid,1,'int');

% Flowtype settings and checks
if ~(fltype < 0 || (fltype>=1 && fltype<=9))
  if pou; fltype=1; end
  if ~pou; fltype=2; end
end

fprintf(['\nStart time: ',num2str(t/dstar),'\n']);

if fltype< -2 || fltype > 9 || fltype==0
  fprintf(['\nInvalid flow type :',num2str(fltype),'\n']);
end

if fltype==1 || fltype==2 || fltype==4 || fltype==5
  dstar = 1;
end

% Move from internal to external quantities
re = re*dstar;
xl = xl/dstar;
zl = zl/dstar;
t = t/dstar;

xs = xs/dstar;

% Write file info
fprintf([' xl= ',num2str(xl),' zl= ',num2str(zl),' t= ',num2str(t),' re= ',num2str(re),'\n'])
if (fltype ==  3); fprintf([' yl= ',num2str(2./dstar),'n\']);        end

fprintf('\nFlow type:');
if (fltype == -2); fprintf('\n Falkner-Skan-Cooke boundary layer\n');   end
if (fltype == -1); fprintf('\n Falkner-Skan boundary layer\n');         end
if (fltype ==  1); fprintf('\n Poiseuille flow\n');                     end
if (fltype ==  2); fprintf('\n Couette flow\n');                        end
if (fltype ==  3); fprintf('\n Blasius boundary layer flow\n');         end
if (fltype ==  4); fprintf('\n spatial Poiseuille flow\n');             end
if (fltype ==  5); fprintf('\n spatial Couette flow\n');                end
if (fltype ==  6); fprintf('\n spatial Blasius boundary layer\n');      end   
if (fltype ==  7); fprintf('\n spatial Falkner-Skan boundary layer\n'); end 
if (fltype ==  8); fprintf('\n spatial Falkner-Skan-Cooke boundary layer\n'); end 
if (fltype ==  9); fprintf('\n spatial parallel boundary layer\n');     end 
 
vvar{1} = 'u';
vvar{2} = 'v';
vvar{3} = 'w';
vvar{4} = 'omx';
vvar{5} = 'omy';
vvar{6} = 'omz';

if (tpl==1); fprintf(['\nxy-plane of ',num2str(vvar{ivar}),' at z= ',num2str(cpl),'\n']); end
if (tpl==2); fprintf(['\nxz-plane of ',num2str(vvar{ivar}),' at y= ',num2str(cpl),'\n']); end
if (tpl==3); fprintf(['\nyz-plane of ',num2str(vvar{ivar}),' at x= ',num2str(cpl),'\n']); end

% File parameters
fprintf(['\nFile parameters nx,nyp,nzc,nfzsym :\n',num2str([nxin,nypin,nzcin,nfzsin]),'\n']);

% This number tells the location in the binary file (i.e. "the byte where")
% the planes start to read.

fstart_plane = ftell(fid);

% Number of bytes per plane, with info about time and xsa (and eol)
% This number can be used to tell the reader to jump times (planes)

if tpl == 1; Nperplane = (64*nxin*nypin + 64*2 +32*4)/8; end
if tpl == 2; Nperplane = (64*nxin*nzcin + 64*2 +32*4)/8; end
if tpl == 3; Nperplane = (64*nypin*nzcin + 64*2 +32*4)/8; end

% Compute the number of planes from the bytes
fseek(fid,0,'eof');
Nplanes = (ftell(fid)-fstart_plane)/Nperplane;

% Actual planes to be read
npl = min(maxNT,Nplanes);

% Print some info
fprintf(['\nNumber of saved planes        : ',num2str(Nplanes),'\n']);
fprintf(['Maximum planes to read (maxNT): ',num2str(maxNT),'\n']);
fprintf(['Actual planes to read         : ',num2str(npl),'\n']);

% Place the pointer back to the first plane available
fclose(fid);

% Initialize containers
if tpl == 1; out = cell(npl,1); for it = 1:npl; out{it} = NaN(nxin,nypin); end; end
if tpl == 2; out = cell(npl,1); for it = 1:npl; out{it} = NaN(nxin,nzcin); end; end
if tpl == 3; out = cell(npl,1); for it = 1:npl; out{it} = NaN(nypin,nzcin); end; end

ta = NaN(npl,1);
xsa = NaN(npl,1);

% Start parallel environment if none exists
if isempty(gcp)
  mycl = parcluster('local');
  % use max number of processes available
  nproc = mycl.NumWorkers;
  parpool(nproc);
end

% I take 10% of the data as reference for printing info
plref = floor(npl/10);
fprintf('\nReading ...\n');

% Start looping over the times
parfor ipl = 1:npl
  if mod(ipl,plref)==0; fprintf('... [~ +10%% done]\n'); end;

  fid = fopen(filename,'r',file_format);
  
  idline = (ipl-1)*Nperplane + fstart_plane;
  fseek(fid,idline,'bof');
  eol = fread(fid,1,'int');
  ta(ipl)  = fread(fid,1,'float64');
  xsa(ipl) = fread(fid,1,'float64');
  eol = fread(fid,1,'int');
  
  eol = fread(fid,1,'int');
  if tpl == 1;
    out{ipl}(:,:) = fftshift(reshape(fread(fid,nxin*nypin,'float64'),[nxin,nypin]),1);
  end
  if tpl == 2;
    out{ipl}(:,:) = fftshift(reshape(fread(fid,nxin*nzcin,'float64'),[nxin,nzcin]),1);
  end
  if tpl == 3;
    out{ipl}(:,:) = reshape(fread(fid,nypin*nzcin,'float64'),[nypin,nzcin]);
  end
  eol = fread(fid,1,'int');
  fclose(fid);
end

% From internal to external quantities
param.ta = ta/dstar;

% Build axes
if tpl == 1
  dx = xl/nxin; param.argx = 0:dx:xl-dx;
  y = chebdif(nypin,1); param.argy = (y(end:-1:1)+1)/dstar;
elseif tpl == 2
  dx = xl/nxin;  param.argx = 0:dx:xl-dx;
  dz = zl/nzcin; param.argy = -zl/2:dz:zl/2-dz;
elseif tpl == 3
  y = chebdif(nypin,1); param.argy = (y(end:-1:1)+1)/dstar;
  dz = zl/nzcin; param.argx = -zl/2:dz:zl/2-dz;
end

param.re = re;
param.t = t;
fprintf('\nDone.\n');


toc
