function [ output_args ] = visual( x,y )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualization
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

if nargin == 1
	x = x(:);
	X = stft(x);
	X = abs(squeeze(X)');
	figure,
	imagesc(10*log10(X))
	set(gca,'YDir','normal')
	caxis([-50 10])
	colorbar
	title('input signal')
	xlabel('frame')
	ylabel('frequeny')
end

if nargin == 2
	X = stft(x);
	Y = stft(y);
	X = abs(squeeze(X)');
	Y = abs(squeeze(Y)');
	figure,
	subplot(211)
	imagesc(10*log10(X))
	set(gca,'YDir','normal')
	caxis([-50 10])
	colorbar
	title('input signal')
	xlabel('frame')
	ylabel('frequeny')
	subplot(212)
	imagesc(10*log10(Y))
	set(gca,'YDir','normal')
	caxis([-50 10])
	colorbar
	title('output signal')
	xlabel('frame')
	ylabel('frequeny')
end



end

