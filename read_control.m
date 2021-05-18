function [output] = read_control(file,content)
%  
%  read_control: reads the control.o file from SIMSON (control plugin)
%  file:         the name of file (usually contorl.o but could be different)
%  content:      the array with the amount of signals saved
%       ASCII:   [numy+numz numu numdio numgio]
%       binary:  [numt numy+numz numu numdio numgio]
%  output:       structured object containing the time series and the signals
%                specified in content
%
%  Pierluigi Morra (pmorra@mech.kth.se)

dot_index = find(file=='.',1,'last');
if length(file(dot_index:end)) == 2
  if file(dot_index:end) == '.o'
    file_type = 'ASCII';
    disp('')
    disp([file,' is a ASCII file.'])
    disp('')
    if length(content) ~= 4
      disp('For ASCII file: "content" length must be 4')
      return
    end
  end
else
  file_type = 'binary';
  disp('')
  disp([file, ' is a binary file.'])
  disp('')
  if length(content) ~= 5
      disp('For binary file: "content" length must be 5')
      return
  end
end


switch file_type
  
  case 'ASCII'
    s = load(file);

    numyy   = content(1); % Number of total measurments (sensors)
    %numz   = content(3); % Number of z-references  (sensors)
    numuu   = content(2); % Number of total inputs      (actuators)
    numdio  = content(3); % Number of input disturbances  (on inputs)
    numgio  = content(4); % Number of output disturbances (on outputs)
    
    
    %if content(1) > 0
    output.t  = s(:,1); last = 2;
    %end
    if numyy > 0
        output.yy = s(:,last:numyy+1);
        last = 1+numyy+1;
    end
%     if numz > 0
%         output.z = s(:,last:numz+last-1);
%         last = last+numz;
%     end
    if numuu > 0 
        output.uu = s(:,last:numuu+last-1);
        last = last+numuu;
    end
    if numdio > 0
        output.dio  = s(:,last:numdio+last-1);
        last = last+numdio;
    end
    if numgio > 0
        output.ngio = s(:,last:numgio+last-1);
        last = last+numgio;
    end
    
  case 'binary'  
    output.t = NaN(content(1),1);
    if content(2) ~= 0; output.yy    = NaN(content(1),content(2)); end
    %if content(3) ~= 0; output.z    = NaN(content(1),content(3)); end
    if content(3) ~= 0; output.uu    = NaN(content(1),content(3)); end
    if content(4) ~= 0; output.dio   = NaN(content(1),content(4)); end
    if content(5) ~= 0; output.ngio  = NaN(content(1),content(5)); end
    
    file_format = 'ieee-le.l64';

    % Open the file and check for endianness
%     fid = fopen(file,'r',file_format);
%     eol=fread(fid,1,'int');
%     
% %     % Check file endianness
% %     if eol ~= 40 
% %         fclose(fid);
% %         file_format = 'ieee-be.l64';
% %         fid = fopen(file,'r',file_format);
% %     end
%     fclose(fid)
    disp(' ')
    disp(['The binary file is opened with the string: "',file_format,'"']);
    disp(' ')
    
    fid = fopen(file,'r',file_format);
    for i = 1:content(1)
      eol = fread(fid,1,'int');
      output.t(i) = fread(fid,1,'float64');
      if content(2) ~= 0; output.yy(i,:)    = fread(fid,content(2),'float64'); end
      %if content(3) ~= 0; output.z(i,:)    = fread(fid,36,'float64'); end
      if content(3) ~= 0; output.uu(i,:)    = fread(fid,content(3),'float64');  end
      if content(4) ~= 0; output.dio(i,:)   = fread(fid,content(4),'float64');  end
      if content(5) ~= 0; output.ngio(i,:)  = fread(fid,content(5),'float64');  end
      eol = fread(fid,1,'int');
    end
    fclose(fid);
end
end