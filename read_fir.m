function [H,nt,ny,dt,iy,t] = read_fir(file)
%
%   [H,nt,ny,dt,iy,t] = read_fir(file)
%

f = fopen(file,'r');

dt = fscanf(f,'%e\t',[1 1]);
nt = fscanf(f,'%d\t',[1 2]);
ny = fscanf(f,'%d\t',[1 2]);
iy = ny(1):ny(2);
%iy = iy(1):iy(2);

H = zeros(ny(2)-ny(1)+1,nt(2)-nt(1)+1);
for i = 1:nt(2)-nt(1)+1
    H(:,i) = fscanf(f,'%e\t',ny(2)-ny(1)+1);
end

f = fclose(f);

t = dt*(nt(1)-1:nt(2)-1);