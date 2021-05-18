function [vel,x,y,z,Lvec,maxt,Re,dstar] = read_flowfield(filename)
%%% Read_flowfield: makes use of readdns, fou2phys, mat2vel and gives back
%%%            the flowfield in the physical space.
%%%  INPUT: (filename)
%%% OUTPUT: [vel,x,y,z,Lvec,maxt,Re,dstar]

[vel,xF,yF,zF,Lx,Ly,Lz,maxt,Re,~,dstar]=readdns(filename);
Lvec=[Lx,Ly,Lz];
vel=fou2phys(vel,0,0);

[vel,x,y,z]=mat2vel(vel,xF,yF,zF);
end

