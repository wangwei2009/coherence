function [ y,Fvv,SNR] = process(x,d,varargin)
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
if nargin==3
    method = varargin{1};
else
    method = 3;
end;

X = stft(x,512,512,256);
N_FFT = (size(X,3)-1)*2;
N = size(x,2);
half_bin = size(X,3);

P.Pxii = ones(N,half_bin);
P.Pxij = ones(N,N,half_bin);
P.Fvv = zeros(N,N,half_bin);

G = zeros(1,half_bin);
SNR = zeros(size(X,1),half_bin);

% output spectral
Y = squeeze(X(:,1,:));

alpha = 0.60;
alpha_MSC = 0;

iNumFilts = 40;
%aad_H     = ComputeFilterResponse(iNumFilts, N_FFT);

for frameIndex = 1:size(X,1)
    X_t = squeeze(X(frameIndex,:,:));
    P = update_MSC(X_t,P,alpha,alpha_MSC);
    
    Gmin = 0.1;
%     method = 5;    %3/5
    for k = 1:half_bin
        [G(k),SNR(frameIndex,k)] = getweights(P.Fvv,k,d,Gmin,method);
%         G(k) = sqrt(G(k));
    end
    ad_X_Bar = squeeze(mean(X(frameIndex,:,:))).';
    G(1:16) = sqrt(G(1:16));
%     [aad_X_tilde] = ChannelWeighting(ad_X_Bar,G,aad_H);
    
%     aad_X_tilde(1:8) = ad_X_Bar(1:8);
%     Y(frameIndex,:) = aad_X_tilde;
    Y(frameIndex,:) = Y(frameIndex,:).*G;
%     Y(frameIndex,:) = squeeze(mean(X(frameIndex,:,:))).'.*G;
    
end
y = istft(Y,512,512,256);
scale = 1.5;
y = real(y)*scale;
Fvv = P.Fvv;

end

