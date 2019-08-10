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
d = 0.02;
r = d/2; 

%% generate diffuse noise
[ noise,Pos] = noise_gen_URA(r);

noise = noise'/500;
N = size(noise,2);        %Channels
%% signal simulation

pathname = '../../sound/';

%use a clean speech audio as desired signal
[speech ,fs] = audioread([pathname,'S_29_04.wav']);
[interference ,fs] = audioread([pathname,'S_72_09.wav']);

               
% speech(20000:30000) = 0.01;
[s,h]  = generate_signal(speech,[0,0],r,1,0.2);
interf = generate_signal(interference,[90,0],r,1,0.2);

% signal+interference+diffuse noise
len_min = min(min(size(s,1),size(interf,1)),size(noise,1));
speech = speech(1:len_min);
switch 1
    case 1
        x = s(1:len_min,:)+interf(1:len_min,:)+noise(1:len_min,:);
    case 2
        x = s(1:len_min,:)+interf(1:len_min,:);
    case 3
        x = s(1:len_min,:);
    otherwise
        disp('other value')
end

slice = [1,3];
x = x(:,slice);
SNR = snr(s(1:len_min,slice),interf(1:len_min,slice))
% x = [speech(3:end),speech(1:end-2)]+[interf(1:size(speech(3:end),1))',interf(1:size(speech(3:end),1))'];
% d = 0.0213*2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = size(x,2);
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

x1 = x;
%% estimate noise coherence function
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
[pesq_mos]= pesq_vec(speech, out,fs)
% visual( x(:,1),out );
% util.fig(out, fs);


