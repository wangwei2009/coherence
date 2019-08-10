%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% endfire
% refer to "A Dual-Microphone Speech Enhancement Algorithm
% Based on the Coherence Function"
%
% broadside
% refer to "A coherence-based noise reduction algorithm for binaural
% hearing aids"
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clear all;
addpath(genpath('../../../../lib'));
c = 340; % speed of sound

%%
%% load recorded office noise audio

fs = 16000;

angle = [0,0]/180*pi;
% array spacing
d = 0.0213;
r = d/2; 

slice = [1,3];
[ sig ] = signal_simulation( r,slice );
x = sig.x;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = size(x,2);
%% 
frameLength = 256;
overlap = 128;
inc = frameLength - overlap;
N_FFT = 256;

%% fixed wideband MVDR using pre-defined noise cross-spectral matrix
Fvv2 = ones(N_FFT/2+1,M,M);     % postfilter MSC
P_len = N_FFT/2+1;
% Pxii_pre = ones(N,P_len);
% Pxij_pre = ones((N*N-N)/2,P_len);

x1 = x;
%% estimate noise coherence function
noise = sig.clean_n;
for i = 1:size(noise,2)
    for j = 1:size(noise,2)       
        [sc,F]=mycohere(noise(:,i),noise(:,j),N_FFT,fs,hanning(N_FFT),0.75*N_FFT);
        Fvv2(:,i,j) = sc;%real(sc);
    end
end

Fvv = GenNoiseMSC(M,N_FFT,fs,d);  % superdirective noise CSD

%% process
[ out,Fvv2,SNR] = process(x1,d);

%% evaluate
speech = sig.speech;
[pesq_mos]= pesq_vec(speech, out,fs)
% visual( x(:,1),out );
% util.fig(out, fs);


