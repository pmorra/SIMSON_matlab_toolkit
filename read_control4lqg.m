function [output] = read_control4lqg(file,content)

% ASCII or BINARY
dot_index = find(file=='.',1,'last');
if length(file(dot_index:end)) == 2
  if file(dot_index:end) == '.o'
    file_type = 'ASCII';
    disp('')
    disp([file,' is a ASCII file.'])
    disp('')
  end
else
  file_type = 'binary';
  disp('')
  disp([file, ' is a binary file.'])
  disp('') 
end

% Read content map
nyio = content(1);
nuio = content(2);
ndio = content(3);
ngio = content(4);
if strcmp(file_type,'binary')
  nt = content(5);
end
  

switch file_type
  
  case 'ASCII'
    data = load(file);

    last = 0;

    % time
    output.t  = data(:,last+1); last = last+1;

    % outputs
    if nyio > 0
        output.yy = data(:,last+(1:nyio)); last = last+nyio;
    end

    % inputs
    if nuio > 0 
        output.uu = data(:,last+(1:nuio)); last = last+nuio;
    end

    % input disturbances
    if ndio > 0
        output.dd = data(:,last+(1:ndio)); last = last+ndio;
    end

    % output disturbances
    if ngio > 0
        output.gg = data(:,last+(1:ngio)); last = last+ngio;
    end
    
  case 'binary'
    % initialize with NaN
    output.t = NaN(nt,1);
    if nyio ~= 0; output.yy = NaN(nt,nyio); end
    if nuio ~= 0; output.uu = NaN(nt,nuio); end
    if ndio ~= 0; output.dd = NaN(nt,ndio); end
    if ngio ~= 0; output.gg = NaN(nt,ngio); end
    
    % input file format
    file_format = 'ieee-le.l64';
    
    disp(' ')
    disp(['The binary file is opened with the string: "',file_format,'"']);
    disp(' ')
    
    % read data
    fid = fopen(file,'r',file_format);
    for i = 1:nt
      eol = fread(fid,1,'int');
      output.t(i) = fread(fid,1,'float64');
      if nyio ~= 0; output.yy(i,:) = fread(fid,nyio,'float64'); end
      if nuio ~= 0; output.uu(i,:) = fread(fid,nuio,'float64'); end
      if ndio ~= 0; output.dd(i,:) = fread(fid,ndio,'float64'); end
      if ngio ~= 0; output.gg(i,:) = fread(fid,ngio,'float64'); end
      eol = fread(fid,1,'int');
    end
    fclose(fid);

end