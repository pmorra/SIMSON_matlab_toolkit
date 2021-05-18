function L2plotter(x,y,z,arg,percentage)
%   Produce a 3D isosurface plot in the simson box of x,y,z.
%   The function is intended for lambda2 criterion vortex visualization.
%   "arg" is the lambda2 value on every point of the box.
%   "vel" is used to color the isosurfaces (example: one could use the
%         velocity value)
%   "percentage" is the percentage of the maximum to compute the
%                isosurfaces and must be in [0,1].
%   "arg" and "vel" must be consistent with x,y,z

arg   = permute(arg(:,:,:),[3 1 2]);
%vel   = permute(vel(:,:,:),[3 1 2]);

xlim(gca,[x(1) x(end)+1]);
ylim(gca,[z(1) z(end)]);
zlim(gca,[y(1) y(end)]);

ma    = 1;%max(max(max(-arg)));

p     = patch(isosurface(x,z,y,arg,-percentage*ma));
%cdata = smooth3(vel,'box',3);

p.FaceColor='red';
p.EdgeColor='none';

isonormals(x(1:end),z(1:end),y(1:end),arg,p);
%isocolors(x(1:end),z(1:end),y(1:end),cdata,p);
%colormap('jet');
daspect([1 1 0.5]);
view(3);camlight; lighting gouraud

end
