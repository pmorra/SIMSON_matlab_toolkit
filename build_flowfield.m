function build_flowfield(filename,vel,Re,Lx,Ly,Lz,flength,varargin)
%
% This function builds a velocity field in SIMSON format, which can be used
% as input in SIMSON directly.
% It assumes the presence of a working copy of SIMSON at the address
% "$HOME/codes/simson". In case it does not exist the function can not be
% used.
% It also assumes the presence of the script "getsimson" in the "$HOME/bin"
% folder. If the script is not there, the function will not work.
% Finally, the SIMSON working copy should have the right bls.f file (email
% to pmorra@mech.kth.se for explanations).
% 
% INPUTS:
% 
% filename: the filename of the output field in SIMSON format.
% vel:      the velocity field as structured object.
%           The "vel" input argument is assumed to be a structured object
%           such as the three components of the velocity are vel.u, vel.v,
%           and vel.w with dimensions (nx,ny,nz). This is also the output
%           from the matlab function read_flowfield(args).
% Re:       Reynolds with respect to displacement thickness at inlet.
% Lx:       streamwise box length.
% Ly:       wall-normal box length (height).
% Lz:       spanwise box length.
%
% flenght:  fringe length.
%
% (c) 2017, Stockholm. Pierluigi Morra.
%
%     email: pmorra@mech.kth.se
%

% parse inputs
write_dat_file = 'n';
already_compiled = 'n';
clean_it = 'y';
% NOTE: write something to parse the "varargin"
if nargin == 8
  write_dat_file = varargin{1};
elseif nargin == 9
  write_dat_file = varargin{1};
  already_compiled = varargin{2};
elseif nargin == 10
  write_dat_file = varargin{1};
  already_compiled = varargin{2};
  clean_it = varargin{3};
elseif nargin > 10
  error('ERROR:\nToo many input arguments.')
end

% load external library
addpath(pathdef);

% check filename
pos = strfind(filename,'.u');
if ~isempty(pos) && pos == length(filename)-1
    filename = filename(1:end-2);
end
pos = strfind(filename,'.dat');
if ~isempty(pos) && pos == length(filename)-3
    filename = filename(1:end-4);
end

% retrieve dimensions
[nx,ny,nz] = size(vel.u);

% permute
u = permute(vel.u,[1 3 2]);
v = permute(vel.v,[1 3 2]);
w = permute(vel.w,[1 3 2]);

% shift along x (streamwise or first dimension) direction
dum1 = u(nx/2+1:end,:,:);
dum2 = u(1:nx/2,:,:);
u(nx/2+1:end,:,:) = dum2;
u(1:nx/2,:,:)     = dum1;

dum1 = v(nx/2+1:end,:,:);
dum2 = v(1:nx/2,:,:);
v(nx/2+1:end,:,:) = dum2;
v(1:nx/2,:,:)     = dum1;

dum1 = w(nx/2+1:end,:,:);
dum2 = w(1:nx/2,:,:);
w(nx/2+1:end,:,:) = dum2;
w(1:nx/2,:,:)     = dum1;

clear dum1 dum2

% write .dat file
fid = fopen([filename,'.dat'],'w');

fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',ny);
fprintf(fid,'%d\n',nz);

for j = ny:-1:1
  fprintf(fid,'%23.16e ',u(:,:,j)); fprintf(fid,'\n');
  fprintf(fid,'%23.16e ',v(:,:,j)); fprintf(fid,'\n');
  fprintf(fid,'%23.16e ',w(:,:,j)); fprintf(fid,'\n');
end

fclose(fid);


% dealiasing flags
nfxd = 1; disp(['Dealiasing in x: nfxd = ',num2str(nfxd)])
nfyd = 0; disp(['Dealiasing in y: nfyd = ',num2str(nfyd)])
nfzd = 1; disp(['Dealiasing in z: nfzd = ',num2str(nfzd)])

% set par.f
disp('Looking for par.f ...')
ishere_par = exist('par.f');
if ~exist('par.f')
  disp('par.f not found. Making one...')
  system('cp $HOME/codes/simson/par.f .');
  system(['sed -i ''s/(nx=32,ny=33,nz=32)/(nx=',...
          num2str(nx),',ny=',num2str(ny),',nz=',num2str(nz),')/g'' par.f']);
  system(['sed -i ''s/(nfxd=1,nfyd=0,nfzd=1)/(nfxd=',...
          num2str(nfxd),',nfyd=',num2str(nfyd),',nfzd=',num2str(nfzd),')/g'' par.f']);
  disp('Done.')
else
  disp('par.f found. Using it.')
end

% set fsc.i, compile, execute
disp('Looking for fsc.i ...')
ishere_fsc = exist('fsc.i');
if ~exist('fsc.i')
  disp('fsc.i not found. Making one ...')
  fid = fopen('fsc.i','w');

  fprintf(fid,'0         m:    power law exponential' ); fprintf(fid,'\n');  
  fprintf(fid,'16392     n:    wall normal resolution'); fprintf(fid,'\n'); 
  fprintf(fid,'1.0e-15   eps:  convergence criterion' ); fprintf(fid,'\n');
  fprintf(fid,'200.      ymax: box height'            ); fprintf(fid,'\n');  
  fprintf(fid,'0         pr:   Prandtl number'        ); fprintf(fid,'\n');    
  fprintf(fid,'0         m1:   scalar exponent'       ); fprintf(fid,'\n');    
    
  fclose(fid);
  disp('Done.')
else
  disp('fsc.i found. Using it.')
end

if already_compiled == 'n'
  system('getsimson "" fsc; ./fsc');
end

% set bls.i, compile, execute
if exist('bls.i') == 2
  system('mv bls.i bls.i.backup');
end
fid = fopen('bls.i','w');

fprintf(fid,[filename,'.u']                                                ); fprintf(fid,'\n');
fprintf(fid,[num2str(Re),'     re:      Reynolds number']                  ); fprintf(fid,'\n'); 
fprintf(fid,[num2str(Lx),'     xlb:     lenght of the box']                ); fprintf(fid,'\n');  
fprintf(fid,[num2str(Ly),'      h2:      height of the box']               ); fprintf(fid,'\n');  
fprintf(fid,[num2str(Lz),'      zlb:     width of the box']                ); fprintf(fid,'\n');     
fprintf(fid,'fsc.dat'                                                      ); fprintf(fid,'\n');     
fprintf(fid,'6        fltype:    base flow type'                           ); fprintf(fid,'\n');     
fprintf(fid,['-',num2str(flength),'    bstart:    fringe starting x']      ); fprintf(fid,'\n');    
fprintf(fid,[num2str(flength),'    bslope:    fringe length']              ); fprintf(fid,'\n');    
fprintf(fid,'.true.   pert:      flag to generate a flow with a base flow' ); fprintf(fid,'\n');    
fprintf(fid,'0.                  tilting angle'                            ); fprintf(fid,'\n');     
fprintf(fid,'0.       ushift:    Galilei shift velocity'                   ); fprintf(fid,'\n');     
fprintf(fid,'.false.  locdi:     flag to generate a localized disturb'     ); fprintf(fid,'\n');        
fprintf(fid,'.false.  gaussian:  flag to generate a gaussian disturb'      ); fprintf(fid,'\n');       
fprintf(fid,'.false.  waves:     flag to generate a pair of oblique waves' ); fprintf(fid,'\n');       
fprintf(fid,'.false.  os:        flag to use tabulated eigenmodes'         ); fprintf(fid,'\n');     
fprintf(fid,'.false.  specm:     spectral space mode'                      ); fprintf(fid,'\n');    
fprintf(fid,'.true.   pertfromfile: flag to read perturbation from file'   ); fprintf(fid,'\n');       
fprintf(fid,'1            amp:   an additional  multiplicative amplitude'  ); fprintf(fid,'\n');    
fprintf(fid,'.true.       pertphys: read from file (filename on next line)'); fprintf(fid,'\n');  
fprintf(fid,[filename,'.dat']                                              ); fprintf(fid,'\n');        
fprintf(fid,'.false.  noise:     flag to add noise'                        ); fprintf(fid,'\n');          

fclose(fid);

if already_compiled == 'n'
  system('getsimson "" bls; ./bls');
elseif already_compiled == 'y'
  system('./bls');
end

% clean folder from unnecessary files (everything set before)
if clean_it == 'y'
  if ~ishere_par
    system('rm par.f');
  end
  if ~ishere_fsc
    system('rm fsc.i fsc.dat');
  end
  if exist('bls.i.backup') == 2
    system('rm bls.i; mv bls.i.backup bls.i');
  else
    system('rm bls.i');
  end
  system('rm fsc bls');
end
if strcmp(write_dat_file,'n')
    system(['rm ',filename,'.dat']);
end
end
