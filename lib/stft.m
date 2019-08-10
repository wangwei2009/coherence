function Y = stft( x,nfft,frameLength,inc,window,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-Channle Short Time Fourier Transform
% 
%
% example Usage: 
%   Y = stft( x,256,256,128)
%
% Inputs:
%   x              input data,stored in column-wise
%   nfft          fft points
%   frameLength    frame length
%   inc            hop size
%   window         analysis window        
%   varargin
%       
%
% Outputs:
%   Y            STFT matrix,[frameNum,channel,half_bin]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x = x(:); % make a column vector for ease later

if nargin == 1
	nfft = 256;
end
if nargin < 3
   frameLength = nfft;
   inc = frameLength/2;
end
if nargin < 5
    window = sqrt(hann(frameLength));
end


M = size(x,2);            % channel
dataLength = size(x,1);   % data length

half_bin = nfft/2+1;

frameNum = fix((dataLength - (frameLength-inc))/inc); % Number of frame

Y = zeros(frameNum,M,half_bin); % outout STFT matrrix

for frameIndex = 1:frameNum
    xt = x((frameIndex-1)*inc+1:(frameIndex-1)*inc+frameLength,:); % enframe
    xt = xt.*window; % windowing
    Xt = fft(xt,nfft);    % fft
    Y(frameIndex,:,:) = Xt(1:half_bin,:).';
end

end

