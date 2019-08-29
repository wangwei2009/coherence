function x = istft( Y,nfft,frameLength,inc,window,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse Short Time Fourier Transform
% 
%
% example Usage: 
%   x = stft( Y)
% 
%
% Inputs:
%   Y              input time-frequency matrix,stored as [frameNum,half_bin]    
%   nfft           fft points
%   frameLength    frame length
%   inc            hop size
%   window         analysis window   
%   varargin
%       
%
% Outputs:
%   x            time-domain data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = squeeze(Y);
if nargin == 1
    nfft = (size(Y,2)-1)*2;                % data length
end
if nargin < 3
   frameLength = nfft;
   inc = frameLength/2;
end
if nargin < 5
    window = sqrt(hann(frameLength+1));
    window = window(1:end-1);
end

frameNum = size(Y,1);            % channel

x = zeros(frameNum*inc+(frameLength-inc),1);

for frameIndex = 1:frameNum
    Yt_half = Y(frameIndex,:);
    Yt = [Yt_half,conj(fliplr(Yt_half(2:end-1)))].';
    xt = ifft(Yt,nfft).*window;     % windowing
    % overlap-add
    x((frameIndex-1)*inc+1:(frameIndex-1)*inc+frameLength) = x((frameIndex-1)*inc+1:(frameIndex-1)*inc+frameLength) + xt;
    
end

end

