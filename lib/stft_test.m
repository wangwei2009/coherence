[x,fs] = audioread('S_01_01.wav');

% x = [zeros(128,1);x];
Y = stft(x);

x_rec = istft(Y);

% Y = stftanalysis(x,256,128);
% 
% x_rec = stftsynthesis(Y,256,128);
