function [ time,x,y,z,Lvec,maxt,Re,dstar,n,fL ]=timings(sbpcheck,subpath,kineticEnergy,velocityVector)
%%%% Timings: gives back a structured object "time" which contains
%%%           all the saved time steps found in the current folder.
%%%  INPUT: (sbpcheck,subpath,kineticEnergy,velocityVector)
%%%           sbpcheck: ['y','n'] checks in the given "subpath" and not 
%%%                     only in the current working directory
%%%            subpath: subpath to check, it is needed if sbpcheck = 'y'
%%%      kineticEnergy: ['y','n'] includes the kinetic energy in "time"
%%%                     output
%%%     velocityVector: ['y','n'] includes the velocity vector in "time"
%%%                     output
%%% OUTPUT: [ time,x,y,z,Lvec,maxt,Re,dstar ]

path=[pwd,'/',subpath];
allFiles = dir([path '/*.u']);
fileList={allFiles.name}; clear allFiles;
n=0; [~,s]=size(fileList);
for i=1:s
    if fileList{i}(1)=='t'
        n=n+1;
    end;
end

FL=zeros(n,1);
fL=cell(n,1); n=0;

n = 0;
for i=1:s
    if fileList{i}(1)=='t';
        if fileList{i}(3)=='_'
            pos = 4;
            n = n+1;
            FL(n) = 0.1*str2num(fileList{i}(pos:end-2)) ...
                    + str2num(fileList{i}(2:pos-2));
        else
            pos = 2;
            n = n+1;
            FL(n) = str2num(fileList{i}(pos:end-2));
        end
    end;
end;

FL = sort(FL);

for i=1:length(FL)
    if FL(i)-floor(FL(i)) > 0
        timeshot = num2str(floor(FL(i)));
        fL{i} = ['t',timeshot,'_'];
        dummy = round(10*(FL(i)-floor(FL(i))));
        if dummy-floor(dummy) > 0
            dummy = round(10*(FL(i)-floor(FL(i))));
        end
        fL{i} = [fL{i},num2str(dummy)];
    else
        timeshot = num2str(FL(i));
        fL{i} = ['t',timeshot];
    end
end

for i=1:n;
    if sbpcheck=='y'
        [time.(fL{i}(1:end)),x,y,z,Lvec,maxt,Re,dstar] = read_flowfield([subpath,'/',fL{i},'.u']); 
    else
        [time.(fL{i}(1:end)),x,y,z,Lvec,maxt,Re,dstar] = read_flowfield([fL{i},'.u']); 
    end
    if kineticEnergy=='y'
    time.(fL{i}(1:end)).('kinEn')=0.5*(time.(fL{i}(1:end)).u.^2+...
        time.(fL{i}(1:end)).v.^2+time.(fL{i}(1:end)).w.^2);
    end
    if velocityVector=='n'
         time.(fL{i}(1:end))=rmfield(time.(fL{i}(1:end)),{'u','v','w'});
    end
end
end