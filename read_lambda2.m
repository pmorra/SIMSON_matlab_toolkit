function [lambda2,status] = read_lambda2(fname)
%
% [l2,status] = read_lambda2(fname)
%

% open the file
status = fopen(fname,'r','ieee-le.l64');
if status < 1
    disp('Error in reading the file')
    return
else
    fid = status;
end

% endianness
eol=fread(fid,1,'int');
if eol ~= 12
  fclose(fid);
  disp(' ')
  disp(['Reading ' fname ' on big endian format'])
  fid=fopen(fname,'r','ieee-be.l64');
  eol=fread(fid,1,'int');
else
  disp(' ')
  disp(['Reading ' fname ' on little endian format'])
end

% read file
nx = fread(fid,1,'int');
ny = fread(fid,1,'int');
nz = fread(fid,1,'int');
eol= fread(fid,1,'int');

% read data
lambda2 = zeros(nx,ny,nz);
for k = 1:nz
    for j = 1:ny
        eol= fread(fid,1,'int');
        lambda2(:,j,k)= fread(fid,nx,'float64');
        eol= fread(fid,1,'int');
    end
end
lambda2=fliplr(lambda2);
status = fclose(fid);
end