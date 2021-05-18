function set_my_pic_style(varargin)
%
% The function looks for all visible figures and modify them according
% to a specified set of settings chosen by Pierluigi Morra.
%
% NOTE: if you do not like the new style, you can always make your own
%       function... Or modify this one.
%


% Look for visible figures
all_figs = findobj('Type','figure');
fprintf('\n')
disp('--------- Set My Picture Style -----------')
if isempty(all_figs)
  % Exit if no figures are found
  disp('Figures found : 0')
  disp('No visible Figures were found! Done.')
  disp('------------------------------------------')
  return
else
  disp(['Figures found : ',num2str(length(all_figs))])
end

my_FS = 20;
if nargin >= 1
  my_FS = varargin{1};
end

% Loop over the figures
for ifig = 1:length(all_figs)
fprintf(['Working on Figure : ',num2str(ifig),' ... '])
h = all_figs(ifig);
for i = 1:length(h.Children)
  % Modify colorbar
  if strcmp(h.Children(i).Type,'colorbar')
    h.Children(i).TickLabelInterpreter = 'latex';
    h.Children(i).FontSize = my_FS;
  end
  % Modify legend
  if strcmp(h.Children(i).Type,'legend')
    lgd = h.Children(i).String;
    for j = 1:length(lgd)
      if lgd{j}(1)   ~= '$'; lgd{j} = ['$' lgd{j}];  end
      if lgd{j}(end) ~= '$'; lgd{j} = [lgd{j} '$'];  end
    end
    h.Children(i).String = lgd;
    h.Children(i).Interpreter = 'latex';
    h.Children(i).FontSize = my_FS;
  end
  % Modify axes
  if strcmp(h.Children(i).Type,'axes')
    h.Children(i).FontSize = my_FS;
    if ~isempty(h.Children(i).XLabel.String)
    if h.Children(i).XLabel.String(1)   ~= '$';  h.Children(i).XLabel.String = ['$' h.Children(i).XLabel.String];  end
    if h.Children(i).XLabel.String(end) ~= '$';  h.Children(i).XLabel.String = [h.Children(i).XLabel.String '$'];  end
    end
    if ~isempty(h.Children(i).YLabel.String)
    if h.Children(i).YLabel.String(1)   ~= '$';  h.Children(i).YLabel.String = ['$' h.Children(i).YLabel.String];  end
    if h.Children(i).YLabel.String(end) ~= '$';  h.Children(i).YLabel.String = [h.Children(i).YLabel.String '$'];  end
    end
    if ~isempty(h.Children(i).ZLabel.String)
    if h.Children(i).ZLabel.String(1)   ~= '$';  h.Children(i).ZLabel.String = ['$' h.Children(i).ZLabel.String];  end
    if h.Children(i).ZLabel.String(end) ~= '$';  h.Children(i).ZLabel.String = [h.Children(i).ZLabel.String '$'];  end
    end
    h.Children(i).Title.Visible = 'off';
    if nargin >= 2 && ~isempty(h.Children(i).Title.String)
      if h.Children(i).Title.String(1)   ~= '$';  h.Children(i).Title.String = ['$' h.Children(i).Title.String];  end
      if h.Children(i).Title.String(end) ~= '$';  h.Children(i).Title.String = [h.Children(i).Title.String '$'];  end
      h.Children(i).Title.Interpreter = 'latex';
      h.Children(i).Title.Visible = varargin{2};
    end
    h.Children(i).TickLabelInterpreter = 'latex';
    h.Children(i).XLabel.Interpreter = 'latex';
    h.Children(i).YLabel.Interpreter = 'latex';
    h.Children(i).ZLabel.Interpreter = 'latex';
  end
end
%set(h,'Position',[0 0 750 665],'PaperPositionMode','auto');
fprintf('Done.\n')
end
disp('Finished.')
disp('------------------------------------------')