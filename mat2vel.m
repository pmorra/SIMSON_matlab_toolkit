function [vel,varargout] = mat2vel(mat,varargin)
%
% vel = mat2vel(mat)
%

nx = size(mat,1);
ny = size(mat,3)/3;
nz = size(mat,2);

%2D
if nz == 2
    mat = mat(:,1,:);
    nz = 1;
end
    
vel.u = zeros(nx,ny,nz);
vel.v = zeros(nx,ny,nz);
vel.w = zeros(nx,ny,nz);

vel.u(:,:,:) = permute(mat(:,:,(ny:-1:1)+0*ny),[1 3 2]);
vel.v(:,:,:) = permute(mat(:,:,(ny:-1:1)+1*ny),[1 3 2]);
vel.w(:,:,:) = permute(mat(:,:,(ny:-1:1)+2*ny),[1 3 2]);

%grid
if nargin == 4
    varargout{1} = varargin{1}(1:end-1);
    varargout{2} = varargin{2}(end:-1:1);
    if nz == 1
        varargout{3} = 0;
    else
        varargout{3} = varargin{3}(1:end-1);
    end
end
    