function [ u,v,w,x,y,z,Lvec,Re ] = parallel_read_flowfield(field_dir,times,varargin)
% [ u,v,w,x,y,z,Lvec,Re ] = parallel_read_flowfield(field_dir,times)
%
% INPUTS  fields_dir: source folder to find the fields
%         times: time series for the saved fields
%
% OUTPUTS u,v,w : cell array containing the field time series for each
%                 velocity component
%         x,y,z : box axes
%         Lvec  : box dimensions
%
% parallel_read_flowfield reads a set of SIMSON fields in parallel and
% gives as output a cell array containing all the fields.
% The advantage is the reading in parallel, which speeds up the reading
% process.
%
% Note 1: The fields need to be named "field.#.u"
% Note 2: It assumes the presence of the function "read_flowfield"
% Note 3: Default use the maximum amount of processes set in MatLab 
%
% Pierluigi Morra (pmorra@mech.kth.se)
%

% Initialize input option
nproc = '';         if nargin >= 3;  nproc     = varargin{1};  end
alpha_DFT = false;  if nargin >= 4;  alpha_DFT = varargin{2};  end
beta_DFT  = false;  if nargin >= 5;  beta_DFT  = varargin{3};  end
turn90 = false;     if nargin >= 6;  turn90    = varargin{4};  end
cut_x = false;      if nargin >= 7;  cut_x     = ~isempty(varargin{5});  end
cut_y = false;      if nargin >= 8;  cut_y     = ~isempty(varargin{6});  end
cut_z = false;      if nargin >= 9;  cut_z     = ~isempty(varargin{7});  end
rem_mean = false;   if nargin >= 10; rem_mean  = varargin{8};  end

if nargin >= 7  && cut_x;  x1 = varargin{5}(1);  x2 = varargin{5}(2);  end
if nargin >= 8  && cut_y;  y1 = varargin{6}(1);  y2 = varargin{6}(2);  end
if nargin >= 9  && cut_z;  z1 = varargin{7}(1);  z2 = varargin{7}(2);  end

% Source directory for fields
fname = [field_dir,'/field.',num2str(times(1)),'.u'];

% Initialize dimensions
[~,X,Y,Z,Lvec,~,Re,~] = read_flowfield(fname);

nt = length(times);
if turn90
  x = Z-Z(1);  y = Y;  z = X-(X(end)+X(2)-X(1))/2;
  NX = length(x);  NY = length(y);  NZ = length(z);
else
  x = X;  y = Y;  z = Z;
  NX = length(x);  NY = length(y);  NZ = length(z);
end

if ~cut_x;  x1 = 1;  x2 = length(x);  end;  nx = length(x(x1:x2));  
if ~cut_y;  y1 = 1;  y2 = length(y);  end;  ny = length(y(y1:y2));  
if ~cut_z;  z1 = 1;  z2 = length(z);  end;  nz = length(z(z1:z2));  

x = x(x1:x2);  y = y(y1:y2);  z = z(z1:z2);

% If remove mean flag : true
if rem_mean && ~isempty(varargin{9}) && isstruct(varargin{9});
  my_mean = varargin{9};
elseif rem_mean
  fprintf('\nProblems with the provided mean. STOP\n');
  return
elseif ~rem_mean
  my_mean.u = zeros(NX,NY,NZ); my_mean.v = my_mean.u; my_mean.w = my_mean.u; 
end

% Initialize dummy container for parfor
u = cell(length(times),1);  v = u;  w = u;

% Set number of processes if given
if ~isempty(nproc); if nargin >= 3;  parpool(nproc);  end

% Print input parameters
disp(['Source directory: "',fname,'"']);
disp('Parameters:')
disp(['            DFT along streamwise direction: ',num2str(alpha_DFT)])
disp(['            DFT along spanwise direction:   ',num2str(beta_DFT)])
disp(['            Turn fields 90deg around wall-normal direction: ',num2str(turn90)])
disp(['            Cut streamwise dimension:  ',num2str(cut_x)])
disp(['            Cut wall-normal dimension: ',num2str(cut_y)])
disp(['            Cut spanwise dimension:    ',num2str(cut_x)])
disp(['            Output X: nx = ',num2str(nx),', x1 = ',num2str(x1),' x2 = ',num2str(x2)]);
disp(['            Output Y: ny = ',num2str(ny),', y1 = ',num2str(y1),' y2 = ',num2str(y2)]);
disp(['            Output Z: nz = ',num2str(nz),', z1 = ',num2str(z1),' z2 = ',num2str(z2)]);
disp(['            Remove mean : ',num2str(rem_mean)]);

% Read and manipulate
parfor it = 1:nt
  fname = [field_dir,'/field.',num2str(times(it)),'.u'];
  if mod(it/nt*100,5) < 1e-2; disp(['... ',num2str(it/nt*100),'%']); end; 
  if ~exist(fname,'file'); disp(['FILE NOT FOUND : ',fname]); end 
  vel = read_flowfield(fname,true);
  if rem_mean; vel.u = vel.u-my_mean.u; vel.v = vel.v-my_mean.v; vel.w = vel.w-my_mean.w;  end
  if turn90
    dum_u = zeros(NX,NY,NZ);  dum_v = dum_u;  dum_w = dum_u;
    for iy = 1:NY
      dum_u(:,iy,:) = flipud( fliplr( squeeze( vel.u(:, iy, :)).'));
      dum_v(:,iy,:) = flipud( fliplr( squeeze( vel.v(:, iy, :)).'));
      dum_w(:,iy,:) = flipud( fliplr( squeeze( vel.w(:, iy, :)).'));
    end
    if alpha_DFT && ~beta_DFT 
      u{it} = fft( dum_w,[],1 ); dum_w = u{it};  u{it} = dum_w( x1:x2, y1:y2, z1:z2 ); dum_w = [];
      v{it} = fft( dum_v,[],1 ); dum_v = v{it};  v{it} = dum_v( x1:x2, y1:y2, z1:z2 ); dum_v = [];
      w{it} = fft(-dum_u,[],1 ); dum_u = w{it};  w{it} = dum_u( x1:x2, y1:y2, z1:z2 ); dum_u = [];
    elseif ~alpha_DFT && beta_DFT
      u{it} = fft( dum_w,[],3 ); dum_w = u{it};  u{it} = dum_w( x1:x2, y1:y2, z1:z2 ); dum_w = [];
      v{it} = fft( dum_v,[],3 ); dum_v = v{it};  v{it} = dum_v( x1:x2, y1:y2, z1:z2 ); dum_v = [];
      w{it} = fft(-dum_u,[],3 ); dum_u = w{it};  w{it} = dum_u( x1:x2, y1:y2, z1:z2 ); dum_u = [];
    elseif ~alpha_DFT && ~beta_DFT
      u{it} =  dum_w( x1:x2, y1:y2, z1:z2 ); dum_w = []; 
      v{it} =  dum_v( x1:x2, y1:y2, z1:z2 ); dum_v = [];
      w{it} = -dum_u( x1:x2, y1:y2, z1:z2 ); dum_u = [];
    elseif alpha_DFT && beta_DFT
      u{it} = fft( fft( dum_w,[],1 ) ,[],3 ); dum_w = u{it};  u{it} = dum_w( x1:x2, y1:y2, z1:z2 ); dum_w = [];
      v{it} = fft( fft( dum_v,[],1 ) ,[],3 ); dum_v = v{it};  v{it} = dum_v( x1:x2, y1:y2, z1:z2 ); dum_v = [];
      w{it} = fft( fft(-dum_u,[],1 ) ,[],3 ); dum_u = w{it};  w{it} = dum_u( x1:x2, y1:y2, z1:z2 ); dum_u = [];
    end
  else
    if alpha_DFT && ~beta_DFT
      u{it} = fft(vel.u,[],1);  u{it} = u{it}( x1:x2, y1:y2, z1:z2 );
      v{it} = fft(vel.v,[],1);  v{it} = v{it}( x1:x2, y1:y2, z1:z2 );
      w{it} = fft(vel.w,[],1);  w{it} = w{it}( x1:x2, y1:y2, z1:z2 );
    elseif ~alpha_DFT && beta_DFT
      u{it} = fft(vel.u,[],3);  u{it} = u{it}( x1:x2, y1:y2, z1:z2 );
      v{it} = fft(vel.v,[],3);  v{it} = v{it}( x1:x2, y1:y2, z1:z2 );
      w{it} = fft(vel.w,[],3);  w{it} = w{it}( x1:x2, y1:y2, z1:z2 );
    elseif ~alpha_DFT && ~beta_DFT
      u{it} = vel.u(x1:x2, y1:y2, z1:z2);
      v{it} = vel.v(x1:x2, y1:y2, z1:z2);
      w{it} = vel.w(x1:x2, y1:y2, z1:z2);
    elseif alpha_DFT && beta_DFT
      u{it} = fft( fft(vel.u,[],1 ),[],3 );  %u{it} = u{it}( x1:x2, y1:y2, z1:z2 );
      v{it} = fft( fft(vel.v,[],1 ),[],3 );  %v{it} = v{it}( x1:x2, y1:y2, z1:z2 );
      w{it} = fft( fft(vel.w,[],1 ),[],3 );  %w{it} = w{it}( x1:x2, y1:y2, z1:z2 );
    end
  end
end
delete(gcp('nocreate'))

end

