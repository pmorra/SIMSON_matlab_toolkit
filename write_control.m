function write_control(file,data)
%
%  write_control(file,data)
%
%  write_control: writes the control.o[.bin] file from SIMSON (control plugin)
%  file:          the name of file (usually control.o but could be different)
%  data:          the array with the amount of signals saved
%
%  Pierluigi Morra (pmorra@mech.kth.se)

dot_index = find(file == '.',1,'last');
if length(file(dot_index:end)) == 2
  if file(dot_index:end) == '.o'
    file_type = 'ASCII';
    disp(' ')
    disp([file,' will be a ASCII file.'])
    disp(' ')
  else
    disp(' ERROR: unknown file extension.')
    return 
  end
elseif length(file(dot_index:end)) == 4
  if file(dot_index:end) == '.bin'
    file_type = 'binary';
    disp(' ')
    disp([file, ' will be a binary file.'])
    disp(' ')
  else
    disp(' ERROR: unknown file extension.')
    return
  end
end


switch file_type
  case 'ASCII'
    fid = fopen(file,'w');
    for time = 1:size(data,1)
      fprintf(fid,'%.16E\t',data(time,:));
      fprintf(fid,'\n');
    end
    fclose(fid);
    
  case 'binary'
    file_format = 'ieee-le.l64';
    disp(['Writing binary format: "',file_format,'"']);
    fid = fopen(file,'w',file_format);
    eol = 8*size(data,2); % could not work in FORTRAN reader !!!
    for time = 1:size(data,1)
      fwrite(fid,eol,'int');
      fwrite(fid,data(time,:),'float64');
      fwrite(fid,eol,'int');
    end
    fclose(fid);
end
end