%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refer to "A Dual-Microphone Speech Enhancement Algorithm
% Based on the Coherence Function"
% endfire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
clear all;
addpath(genpath('../../../../lib'));
c = 340; % speed of sound

%%
M = 2;
N = M;
fs = 16000;

angle = [0,0]/180*pi;

%% 
frameLength = 256;
overlap = 128;
inc = frameLength - overlap;
N_FFT = 256;

%% fixed wideband MVDR using pre-defined noise cross-spectral matrix
Fvv2 = ones(N_FFT/2+1,M,M);     % postfilter MSC
P_len = N_FFT/2+1;
Pxii_pre = ones(N,P_len);
Pxij_pre = ones((N*N-N)/2,P_len);

if(0)
    r = 0.0125*2;
    AudioPath = '../wav/631/1/';
    x1 = audioread([AudioPath,'1.wav'])*10;
    x2 = audioread([AudioPath,'2.wav'])*10;
else
    r = 0.032*2;
    AudioPath = '../wav/xmos/rec1/“ÙπÏ-';
    x1 = audioread([AudioPath,'4.wav']);
    x2 = audioread([AudioPath,'2.wav']);
end

x = [x1,x2];
x1 = x;
%% estimate noise coherence function
noise = x1(1:12000,:);
% noise = x1(15000:30000,:);
% noise = x(36000:end,:);
% x = x(40000:end,:);
for i = 1:size(noise,2)
    for j = 1:size(noise,2)
        
        [sc,F]=mycohere(noise(:,i),noise(:,j),N_FFT,fs,hanning(N_FFT),0.75*N_FFT);
        Fvv2(:,i,j) = sc;%real(sc);

    end
end
Fvv2(1,:,:) = 0.90;   

Fvv = GenNoiseMSC(M,N_FFT,fs,r,'linear');  % superdirective noise CSD
% Fvv2 = ones(N_FFT/2+1,M,M)*0.998;     % postfilter MSC
%% fixed wideband MVDR using pre-defined noise cross-spectral matrix
% [ MVDR_out,H,DI,WNG] = superdirectiveMVDR(x,fs,N_FFT,N_FFT,N_FFT/2,r*2,angle);
angle = [90,0]/180*pi;
[ out,Fvv2,SNR] = process(x1,r);
% SNR = fliplr((SNR)');
plotSNR( SNR )
visual( x(:,1),out );
% Cxy = mscohere(noise(:,1),noise(:,2),hanning(N_FFT),0.75*N_FFT,N_FFT,fs);


