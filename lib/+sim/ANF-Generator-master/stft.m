function Y=stft(x,nfft,dM,dN,wintype)
% stft : Short Time Fourier Transform
% ***************************************************************@
% Inputs: 
%    x,     	signal;
%    nfft,  	window length;
%    dM,			sampling step in Time;
%    dN,			sampling step in Frequency;
%    wintype,	window type;
% Usage:
%    Y=stft(x,nfft,dM,dN,wintype);
% Defaults:
%    wintype='Hanning';
%    dN = 1;
%    dM = 0.5*nfft;
%    nfft = minimum of 256 and the length of the signal;
% Notes:
%    B = STFT(A,NFFT,dM,dN,WINTYPE) calculates the STFT
%    for the signal in vector A.  The signal is split into
%    overlapping segments, each of which are windowed and then
%    Fourier transformed to produce an estimate of the
%    short-term frequency content of the signal.

% Copyright (c) 2000. Dr Israel Cohen.
% All rights reserved. Created  17/12/00.
% ***************************************************************@

x = x(:); % make a column vector for ease later
nx = length(x);
if nargin == 1
	nfft = min(nx,256);
end
if nargin < 3
   dM = 0.5*nfft;
   dN = 1;
end
if nargin < 5
	wintype = 'Hamming';
end

if exist(wintype)
   wins=eval([lower(wintype),sprintf('(%g)',nfft)]);
else
   error(['Undefined window type: ',wintype])
end

% find analysis window for the above synthesis window
win=biorwin(wins,dM,dN);

% make an artificial delay at the beginning an end of the input signal x[n]
% in order to make an accurate STFT-ISTFT operation
delay1 = nfft-1;
delay2 = nfft-1;
x = [zeros(delay1,1); x; zeros(delay2,1)];
nx = length(x);

% figure out number of columns for offsetting the signal
% this may truncate the last portion of the signal since we'd
% rather not append zeros unnecessarily - also makes the fancy
% indexing that follows more difficult.
 ncol = fix((nx-nfft)/dM+1);
y = zeros(nfft,ncol);

% now stuff x into columns of y with the proper offset
% should be able to do this with fancy indexing!
colindex = 1 + (0:(ncol-1))*dM;
rowindex = (1:nfft)';
y(:) = x(rowindex(:,ones(1,ncol))+colindex(ones(nfft,1),:)-1);

% apply the window to the array of offset signal segments
y(:) = win(:,ones(1,ncol)).*y;

% now fft y which does the columns
N = nfft/dN;
for k=1:dN-1
   y(1:N,:) = y(1:N,:) + y((1:N)+k*N,:);
end
y = y(1:N,:);
y = fft(y);
if ~any(any(imag(x)))
	y = y(1:N/2+1,:);
end
% take abs, and use image to display results
if nargout == 0
	imagesc(abs(y));axis xy
else
	Y = y;
end