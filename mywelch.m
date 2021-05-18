function [ wmsfft,A ] = mywelch( signal,tspan,n_bins,varargin )
% mywelch: perform welch method FFT
%   NOTE: "in" and "tspan" must be consistent to each others
%         Default Npad = 1024;
%start = length(t)-5*3^2*2^8+1;
%tspan = t(start:end);
nt    = length(tspan);
%in    = unc.y(start:end,posSens);

nsnaps = nt;
%n_bins = 8;
no_bins = 2*n_bins-1;
Npad = 1024;
if nargin == 4
    Npad = varargin{1};
elseif nargin > 4
    fprintf('\nToo many input arguments.\n');
end
    


if floor(nsnaps/n_bins) ~= nsnaps/n_bins
  disp('nsnaps needs to be dividable by no_bins')
  return
end
 
NN = nsnaps/n_bins;
if floor(NN/2)~=NN/2
  disp('NN needs to be even')
  return
end


wmsfft = zeros(Npad,size(signal,2)); %Welch Method signal FFT
A = zeros(Npad,size(signal,2)); %Welch Method signal FFT Absolute value

for i = 1:size(signal,2)
    for j=1:no_bins
      g = signal(1+(j-1)*NN/2:NN+(j-1)*NN/2,i);

      % Hanning
      h=hann(NN);
      norm=mean(h.^2);
      g = (g.*h)/sqrt(norm);
      wmsffthat = fft(g,Npad)/Npad;   
      Ahat = abs(wmsffthat).^2;

      % Add
      A(:,i) = A(:,i) +Ahat/no_bins;
      wmsfft(:,i) = wmsfft(:,i)+wmsffthat/no_bins;
    end
end
end

