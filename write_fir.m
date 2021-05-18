function write_fir(file,H,t,varargin)
% write_fir(filename,H,t,varargin)

if isempty(varargin)
    ny = size(H,1); ny = [-1 1]*floor((ny-1)/2);
else
    ny = varargin{1};
end

dt = t(2)-t(1); nt = round(t([1,end])/dt)+1;

f = fopen(file,'w');

fprintf(f,'%e\t',dt); fprintf(f,'\n');
fprintf(f,'%d\t',nt); fprintf(f,'\n');
fprintf(f,'%d\t',ny); fprintf(f,'\n');

for i = 1:length(t)
    fprintf(f,'%e\t',H(:,i));
    fprintf(f,'\n');
end

f = fclose(f);