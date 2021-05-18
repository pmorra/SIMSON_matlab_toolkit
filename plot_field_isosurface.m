function plot_field_isosurface(x,y,z,arg,Re,percentage,color,varargin)
% A 3D plot of the "arg" value in the simson box of x,y,z.
%   arg must be consistent with x,y,z
% inputs: (x,y,z,arg,Re,percentage,color)

%[Rex,Re_delta,x,delta,x_0]=calc_rex(x,Re);

arg = permute(arg(:,:,:),[3 1 2]);
arg = flip(arg,1);

xlim(gca,[x(1) x(end)]);
ylim(gca,[z(1) ceil(z(end))]);
zlim(gca,[y(1) y(end)]);

[Rex,Re_delta,xa,delta,x_0] = calc_rex(x,Re);

posit = [0,0,2000,1000];
set(gcf,'Position',posit);

p=patch(isosurface(x,z,y,arg,percentage));
p.FaceColor=color ;p.EdgeColor='none';
if nargin == 8
  p.SpecularStrength = varargin{1};
end
daspect([1 1 0.5]);
isonormals(x(1:end),z(1:end),y(1:end),arg,p); %axis image
view(3); set(gca,'Ytick',[z(1) 0 ceil(z(end))],'Ydir','reverse');
camlight; lighting gouraud;

box on;
xlabel('x');ylabel('z');zlabel('y');

end

