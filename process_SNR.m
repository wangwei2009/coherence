function [ y,Fvv,SNR,SNR_predict,P_s,P_i,P_x,Fvv_all_est] = process_SNR(s,interf,d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency domain processing function
%   
%      input :
%          x : input signal ,samples * channel
%          fs: sample rate
%          N : fft length,frequency bin number
%frameLength : frame length,usually same as N
%        inc : step increment
%          d : array element spacing
%      angle : incident angle
%
%     output :
%    MVDRout : MVDR output
%         x1 : presteered signal,same size as x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_source = stft(s);
X_interf = stft(interf);
X_noisy = stft(s+interf);
N = size(s,2);
half_bin = size(X_source,3);

frames = size(X_source,1);

P_s.Pxii = ones(N,half_bin);
P_s.Pxij = ones(N,N,half_bin);
P_s.Fvv = zeros(N,N,half_bin);
P_s.Fvv_all = zeros(frames,N,N,half_bin);

P_i.Pxii = ones(N,half_bin);
P_i.Pxij = ones(N,N,half_bin);
P_i.Fvv = zeros(N,N,half_bin);
P_i.Fvv_all = zeros(frames,N,N,half_bin);

P_x.Pxii = ones(N,half_bin);
P_x.Pxij = ones(N,N,half_bin);
P_x.Fvv = zeros(N,N,half_bin);
P_x.Fvv_all = zeros(frames,N,N,half_bin);

Fvv_all_est = zeros(frames,N,N,half_bin);

G = zeros(1,half_bin);

% output spectral
Y = squeeze(X_noisy(:,1,:));

SNR_L = zeros(size(X_source,1),half_bin);
SNR_R = zeros(size(X_source,1),half_bin);
SNR = zeros(size(X_source,1),half_bin);
SNR_predict = zeros(size(X_source,1),half_bin);

fs = 16000;
c = 340;
f = 0:fs/256:fs/2;

N_FFT = 256;
k = 0:1:N_FFT/2;
% k = k(2:end);
omega = 2*pi*k/N_FFT;

% omega = 2*pi*f/fs;
tau = fs*d/c;

alpha = 0.8;
alpha_MSC = 0;
% X_interf = X_source;
for frameIndex = 1:size(X_source,1)
    
    SNR_L(frameIndex,:) = real(P_s.Pxii(1,:)./P_i.Pxii(1,:));
    SNR_R(frameIndex,:) = real(P_s.Pxii(2,:)./P_i.Pxii(2,:));
    SNR(frameIndex,:) = real(P_s.Pxii(1,:)./P_i.Pxii(1,:));
    
    X_s = squeeze(X_source(frameIndex,:,:));
    P_s = update_MSC(X_s,P_s,alpha,alpha_MSC);
    P_s.Fvv_all(frameIndex,:,:,:) = P_s.Fvv;
    
    X_i = squeeze(X_interf(frameIndex,:,:));
    P_i = update_MSC(X_i,P_i,alpha,alpha_MSC);
    P_i.Fvv_all(frameIndex,:,:,:) = P_i.Fvv;
    
    X_x = squeeze(X_noisy(frameIndex,:,:));
    P_x = update_MSC(X_x,P_x,alpha,alpha_MSC);
    P_x.Fvv_all(frameIndex,:,:,:) = P_x.Fvv;
    
    
%     Fvv_all_est(frameIndex,1,2,:) = squeeze(P_s.Fvv(1,2,:)).'.*SNR(frameIndex,:)./(SNR(frameIndex,:)+1)+...
%                                     squeeze(P_i.Fvv(1,2,:)).'./(SNR(frameIndex,:)+1);
    Fvv_all_est(frameIndex,1,2,:) = exp(1j*omega*tau).*SNR(frameIndex,:)./(SNR(frameIndex,:)+1)+...
                                    1./(SNR(frameIndex,:)+1);
    for k = 1:half_bin
        [G(k),SNR_predict(frameIndex,k)] = getweights4(P_x.Fvv,k,d);
    end
%     G = sqrt(SNR(frameIndex,k)./(1+SNR(frameIndex,k)));
    
    Y(frameIndex,:) = Y(frameIndex,:).*G;
   
end
y = istft(Y);
Fvv = P_s.Fvv;

end

