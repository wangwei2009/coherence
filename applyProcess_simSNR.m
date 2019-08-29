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
clear all;
addpath(genpath('lib'));

% array spacing
d = 0.0213;
r = d/2; 

slice = [1,3];
[ sig ] = sim.signal_simulation( r,slice );

ang = [0,0]*pi/180;

% SNR = snr(clean_s(:,1),clean_i(:,1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = size(sig.x,2);
%% process
frameLength = 256;
overlap = 128;
inc = frameLength - overlap;
N_FFT = 256;

fs = 16000;
c = 340;
f = 0:fs/256:fs/2;

N_FFT = 256;
k = 0:1:N_FFT/2;
omega = 2*pi*k/N_FFT;

tau = fs*d/c;
Fx1x2 = exp(1j*omega*tau);

[ out,Fvv2,SNR,SNR_predict,P_s,P_i,P_x,Fvv_all_est] = process_SNR(sig.clean_s,sig.clean_i,d);
%% compare predict SNR with true SNR
figure,
subplot(2,1,1),plotSNR( SNR ),title('true SNR'),
subplot(2,1,2),plotSNR( SNR_predict ),title('predict SNR')
frame = 280;
figure,plot(10*log10(SNR(frame,:)))
hold on ,plot(10*log10(SNR_predict(frame,:))),legend(['true SNR at frame ',int2str(frame)],...
                                                     ['predict SNR at frame ',int2str(frame)])

Fvv_s_n = squeeze(P_s.Fvv_all(frame,1,2,:));
figure,
plot(real(Fvv_s_n))
hold on,plot(cos(omega*tau*cos(ang(1)))),title('real part of coherence'),
legend('estimated','theoretical');


%% compare predict coherence with true coherence
k = 16;
Fvv_s_k = squeeze(P_s.Fvv_all(:,1,2,k));
Fvv_i_k = squeeze(P_i.Fvv_all(:,1,2,k));
Fvv_x_k = squeeze(P_x.Fvv_all(:,1,2,k));
Fvv_all_est_k = squeeze(Fvv_all_est(:,1,2,k));
figure,
subplot(2,1,1),plot(real(Fvv_x_k)),hold on,plot(real(Fvv_all_est_k))  % magnitude
legend(['coherence of noisy signal at  ',int2str(fix(k*fs/N_FFT)),'Hz'],...
       ['coherence of model at ',int2str(fix(k*fs/N_FFT)),'Hz']),title('magnitude');
subplot(2,1,2),plot(angle(Fvv_x_k)),hold on,plot(angle(Fvv_all_est_k)),title('phase'); % phase

%% estimate coherence affected by signal or noise with SNR
Fvv_s_n = squeeze(P_s.Fvv_all(frame,1,2,:));
Fvv_i_n = squeeze(P_i.Fvv_all(frame,1,2,:));
Fvv_x_n = squeeze(P_x.Fvv_all(frame,1,2,:));
figure,plot(real(Fvv_s_n)),hold on,plot(real(Fvv_i_n)),hold on,plot(real(Fvv_x_n));
hold on,plot(10*log10(SNR(frame,:))/max(abs(10*log10(SNR(frame,:)))))
legend(['coherence of source at frame  ',int2str(frame)],...
       ['coherence of noise at frame ',int2str(frame)],...
       ['coherence of noisy at frame ',int2str(frame)],...
       ['normalized true SNR at frame ',int2str(frame)])
figure,plot(10*log10(SNR(frame,:)));title(['true SNR(dB) at frame ',int2str(frame)])

speech = sig.speech;
% [pesq_mos]= pesq_vec(speech, out,fs)
% visual( x(:,1),out );
rmpath(genpath('lib'));

