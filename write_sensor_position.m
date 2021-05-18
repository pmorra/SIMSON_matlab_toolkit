function write_sensor_position(filename,Lz,nz,ns,x0,varargin)
% Lz = 50;
% nz = 108;
dz = Lz/nz;
z  = -25:dz:Lz/2-dz;

% ns = 36;
dzs = Lz/ns;

errEXP = 5;
if nargin > 5
    errEXP = varargin{1}
end
sprintf('Computing zscale value based on errEXP = %d ...',errEXP)
zs = dzs/2/sqrt(log(10^errEXP));
sprintf('zscale = %0.5f',zs)

z0 = -Lz/2:dzs:Lz/2-dzs;

for i = 1:length(z0)
  z1 = z-z0(i);
  plot(z,exp(-(z1/zs).^2)); hold on
end
xlabel('z'); ylabel('exp(-(z/zs)^2))');

%% write out
% filename = 'list_sensor.dat';
%x0 = 100;
f = fopen(filename,'w');
fprintf(f,'zscale = %0.5f\n',zs);
for i = 1:length(z0)
    H = [x0, z0(i)];
    fprintf(f,'%e\t',H);
    fprintf(f,'  xsh  zsh');
    fprintf(f,'\n');
end
f = fclose(f);