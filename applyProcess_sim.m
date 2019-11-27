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
addpath(genpath('lib'));
addpath(genpath('E:\work\matlab\ehabets\ANF-Generator-master'));
c = 340; % speed of sound

%%
%% load recorded office noise audio

fs = 16000;

angle = [0,0]/180*pi;
% array spacing
d = 0.025;
r = d/2; 

switch 1
    case 1
        slice = [1,3]; % extract speaker-1
        disp('speaker-1 is in front of mic1')
    case 2
        slice = [2,4]; % extract speaker-2
        disp('speaker-2 is in front of mic2')
    otherwise
        disp('other value')
end

[ sig ] = sim.signal_simulation( r,slice );
rmpath(genpath('E:\work\matlab\ehabets\ANF-Generator-master'));
x = sig.x;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = size(x,2);
%% 
frameLength = 256;
overlap = 128;
inc = frameLength - overlap;
N_FFT = 256;

x1 = x;

%% process

[ y,Fvv2,SNR] = process(x1,d,7);

%% evaluate
speech = sig.speech;
% [pesq_mos]= pesq_vec(speech, out,fs)
rmpath(genpath('lib'));
stoi(sig.clean_i(:,1),x(:,1),fs)        %% STOI for noisy speech
stoi(sig.clean_i(1:length(y),1),y,fs)   %% STOI for processed speech
visual( x(:,1),y );
% util.fig(out, fs);


