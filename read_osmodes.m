function [eigfun1,eigfun2,eigfun3,eigfun4,alfa,beta,gamma,omega,ymax] = read_osmodes(filename)
%    =========================================================
%    Read the Orr-Sommerfield modes file needed for
%    Free-Stream Turbulence in SIMSON.
%    (Brandt L., Schlatter P., Henningson D.S. (2004) can be
%     checked for reference)
% 
%    Outputs the eigenfunctions and the modes.
%
%    This version assumes the binary "filename" is in
%    little-endian format and should be able to check
%    for the endianness.
%
%    (c) Pierluigi Morra, Ardeshir Hanifi
%        Linne` FLOW Centre, KTH Royal Institute of Technology
%        Stockholm, May 2017
%    =========================================================

% The default assumption is a little endian file format
file_format = 'ieee-le.l64';

% Open the file and check for endianness
fid = fopen(filename,'r',file_format);

eol=fread(fid,1,'int');

% Check file endianness
if eol ~= 40 
    fclose(fid);
    file_format = 'ieee-be.l64';
    fid = fopen(filename,'r',file_format);
    eol=fread(fid,1,'int');
end
disp(' ')
disp(['The binary file is opened with the string: "',file_format,'"']);
disp(' ')

% Read the file a first time to get some loop parameters
time = fread(fid,40,'uint8=>char');
eol=fread(fid,2,'int');

vers = fread(fid,1,'int'); eol=fread(fid,2,'int');
n    = fread(fid,1,'int')
re   = fread(fid,1,'float64')
ymax = fread(fid,1,'float64')
eol=fread(fid,2,'int');

kkx    = fread(fid,1,'float64');
alfa   = fread(fid,1,'float64') + 1i*fread(fid,1,'float64');
beta   = fread(fid,1,'float64');
gamma  = fread(fid,1,'float64') + 1i*fread(fid,1,'float64');
eol=fread(fid,2,'int');

z1     = fread(fid,1,'int')

fclose(fid);

% Initialize the arrays
eigfun1 = zeros(n+1,z1);
eigfun2 = zeros(n+1,z1);
eigfun3 = zeros(n+1,z1);
eigfun4 = zeros(n+1,z1);

% Read the file a second time to read everything needed.
fid = fopen(filename,'r',file_format);
eol=fread(fid,1,'int');
time=fread(fid,40,'uint8=>char');
eol=fread(fid,1,'int');

 for j = 1:z1
     eol=fread(fid,1,'int');
     % Line j+1
     vers = fread(fid,1,'int');
     eol=fread(fid,2,'int');
     % Line j+2
     n    = fread(fid,1,'int');
     re   = fread(fid,1,'float64');
     ymax = fread(fid,1,'float64');
     eol=fread(fid,2,'int');
     % Line j+3
     kkx(j)    = fread(fid,1,'float64');
     alfa(j)   = fread(fid,1,'float64') + 1i*fread(fid,1,'float64');
     beta(j)   = fread(fid,1,'float64');
     gamma(j)  = fread(fid,1,'float64') + 1i*fread(fid,1,'float64');
     eol=fread(fid,2,'int');
     % Line j+4
     z1    = fread(fid,1,'int');
     modus = fread(fid,1,'int');
     eol=fread(fid,2,'int');
     % Line (j-1)+5
     scale(j) = fread(fid,1,'float64');
     eol=fread(fid,2,'int');
     % Line j+6
     arg = fread(fid,[9,n+1],'float64');
     eol=fread(fid,2,'int');
     % Storing "arg"...
     eta          = arg(1,:)';
     eigfun1(:,j) = arg(2,:)'+1i*arg(3,:)';
     eigfun2(:,j) = arg(4,:)'+1i*arg(5,:)';
     eigfun3(:,j) = arg(6,:)'+1i*arg(7,:)';
     eigfun4(:,j) = arg(8,:)'+1i*arg(9,:)';
     % Line j+7
     eof = fread(fid,1,'float64') + 1i*fread(fid,1,'float64');
     eol=fread(fid,1,'int');   
 end
 fclose(fid);

 omega = kkx;
end
 