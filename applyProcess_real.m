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
close all
% clear all;
%addpath(genpath('lib'));
c = 340; % speed of sound

%%
%% load recorded office noise audio

fs = 16000;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
frameLength = 256;
overlap = 128;
inc = frameLength - overlap;
N_FFT = 256;
% test xmos 4-mic circular array recordings
x = loadwav('wav/xmos/rec/');
d = 0.064;
switch 1
    case 1
        x = x(:,[1,3]); % extract speaker-1
        disp('speaker-1 is in front of mic1')
    case 2
        x = x(:,[4,2]); % extract speaker-2
        disp('speaker-2 is in front of mic4')
    otherwise
        disp('other value')
end

% x = pcmread('../../../pcm/sound/631/test_631_front3/ETAudioDump_000.pcm',2);
% % x = loadpcm('../sound/631/test_659_front/1/');
% x = [x(:,2),x(:,1)]*10;
% x = [x(1:1600000,2),x(1:1600000,1)];
% d = 0.025;
% fs = 16000;
% x = downsample_48k_8k(x);
x = pcmread('x.pcm',2);
x = x(1:800000,:);
x = [zeros(160,2);x];
d = 0.025;
x1 = x;

%% process

[ y,Fvv2,SNR] = process(x1,d,7);

%% evaluate
% speech = sig.speech;
% [pesq_mos]= pesq_vec(speech, out,fs)
%rmpath(genpath('lib'));
visual( x(:,1),y );
% util.fig(out, fs);


