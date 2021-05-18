function varargout = plot_fir(file,label,varargin)
%
%   plot_fir(filename,label[,color])
%

%default values
color = 'k';
lnst  = '-';
ds = 1;
angle = [25,70];
scale = 1;
% options
if nargin > 2
    if nargin == 3
        color = varargin{1};
    else
        for i = 1:2:length(varargin)
            if strcmpi(varargin{i},'Color')
                color = varargin{i+1};
            elseif strcmpi(varargin{i},'DownSample')
                ds = varargin{i+1};
            elseif strcmpi(varargin{i},'LineStyle')
                lnst = varargin{i+1};
            elseif strcmpi(varargin{i},'View')
                angle = varargin{i+1};
            elseif strcmpi(varargin{i},'Rescale')
                scale = varargin{i+1};
            else
                disp(['''',varargin{i},''' is not a valid option.']);
            end
        end
    end
else
    color = 'k';
end

% read FIR data
[P,~,~,~,m,t] = read_fir(file);

% downsampling
it = [1:ds:length(t),length(t)];
t = t(it);
P = P(:,it)/scale;

if length(m) > 1

    % plot each impulse risponse
    Pax = [-1 1]*max([max(abs(P)) eps]);

    for i = 1:size(P,1)
        tpatch = [t      t(end)   t(1)  ];
        Ppatch = [P(i,:) Pax(1)   Pax(1)];
        patch(tpatch,m(i)+0*tpatch,Ppatch,[1 1 1],'EdgeColor','none'), hold on
        plot3(t,m(i)+0*t,P(i,:),lnst,'Color',color), hold on
    end

    % highlight central impulse response (m = 0)
    hh = plot3(t,0*t,P(m==0,:),lnst,'Linewidth',1.5,'Color',color); hold off

    % adjust view
    view(angle); grid on; axis tight

    % axis labels
    xlabel('i  \Deltat'); ylabel('m');
    zlabel([label,'_{m}(i)']); drawnow

else
    
    hh = plot(t,P,lnst,'Linewidth',1.5,'Color',color);
    grid on; axis tight

    % axis labels
    xlabel('i  \Deltat');
    ylabel([label,'(i)']); drawnow
    
end

if nargout > 0
   varargout{1} = hh; 
end
    