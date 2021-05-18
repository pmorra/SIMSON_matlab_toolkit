function write_foulist(Ti,Tf,fs,varargin)
% write_foulist: creates the file to be used in fou for the FFT in time
%                of multiple fields (SIMSON only)
% INPUTS:
%               Ti = initial timestep of sampling
%               Tf = finel timestep of sampling
%               fs = sampling frequency

%Default saved filename
filename = 'list.dat';

%Options
if nargin == 4
    filename = varargin{1};
end

T = Tf-Ti;

f = fopen(filename,'w');
fprintf(f,[num2str(T/fs),'\n']);

for i = 0:fs:T
    fieldname = ['field.',num2str(Ti+i),'.u'];
    fprintf(f,[fieldname,'\n']);
end
fclose(f);