function [content_ctr,content_unc,i_signals_unc,file_unc,file_Pzu,file_K,file_Pyu,n_lines] = read_mimo(file_mimo)
%
% Read parameters from mimo.i file as used in MIMMO software
% See MIMMO documentations for details.
% 
% [content_ctr,content_unc,file_unc,file_Pzu,file_K,file_Pyu,n_lines] = read_mimo(file_mimo)

fid = fopen(file_mimo,'r');
isfound = false;
isend = false;
Pyu_flag = false;
n_lines = 0;
line = [];
content_unc = [];
content_ctr = [];
file_Pyu = [];
file_K = [];
file_Pzu = [];
file_unc = [];
while ~isend

  file_Pyu = file_K;
  file_K = file_Pzu;
  file_Pzu = file_unc;
  file_unc = line;

  n_lines = n_lines +1;
  
  line = fgetl(fid);
  if line == -1
    isend = true;
    if ~Pyu_flag
      file_Pyu = [];
    end
    break
  end
  
  if ~isempty(strfind(line,'content:')) && isempty(strfind(line,'signal:'))
    content_unc = str2num(line(1:strfind(line,'content:')-1));
  elseif ~isempty(strfind(line,'signal:'))
    signals_ctr = str2num(line(1:strfind(line,'signal:')-1));
    nt = signals_ctr(2)-signals_ctr(1)+1;
    ny = signals_ctr(4)-signals_ctr(3)+1;
    nz = signals_ctr(6)-signals_ctr(5)+1;
    nu = signals_ctr(8)-signals_ctr(7)+1;
    nd = signals_ctr(10)-signals_ctr(9)+1;
    ng = signals_ctr(12)-signals_ctr(11)+1;
    content_ctr = [nt ny+nz nu nd ng];
  elseif ~isempty(strfind(line,'yu_contribute:')) && ~isempty(strfind(line,'.true.'))
    if strfind(line,'.true.') == 1
      Pyu_flag = true;
    end
  end
  if ~isempty(content_unc) && ~isempty(content_ctr)
    isfound = true;
  end

end
fclose(fid);

if ~isfound
  disp('  ')
  disp('WARNING: content line for input or output signals was not found.')
  disp('  ')
else
  if ~isempty(content_ctr) && ~isempty(strfind(file_unc(end-4:end),'.o'))
    content_ctr = content_ctr(2:end);
  end
  if ~isempty(content_unc) && ~isempty(strfind(file_unc(end-4:end),'.o'))
    content_unc = content_unc(2:end);
  end
end

% Remove leading/trailing apostrophe
file_Pyu = file_Pyu(2:end-1);
file_K = file_K(2:end-1);
file_Pzu = file_Pzu(2:end-1);
file_unc = file_unc(2:end-1);
i_signals_unc = signals_ctr;

end