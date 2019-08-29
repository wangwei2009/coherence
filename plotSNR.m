function [ output_args ] = plotSNR( SNR )
% plot SNR 
%function
%
% Usage: (1) plotSNR( SNR );  
%
% Inputs:
%   SNR            SNR ,[frameIndex,frequency]
%   N_FFT          fft points
%   fs             sample frequency in Hz
%   r              array aperture
%   varargin
%       arrayType  array type,linear or circular,default linear
%
% Outputs:
%   Fvv            [half_bin,M,M] coherence function
SNR = SNR';
% figure,imagesc([1:size(SNR,1)]*128,[1:size(SNR,2)],10*log10(flipud(abs(SNR))));
imagesc([1:size(SNR,2)]*128,[1:size(SNR,1)]*16000/256,10*log10(abs(SNR)));
set(gca,'YDir','normal')
title('SNR')
xlabel('sample')
ylabel('frequeny')
caxis([-10 10])
colorbar
end

