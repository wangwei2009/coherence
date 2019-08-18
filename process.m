function [ y,Fvv,SNR] = process(x,d)
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

X = stft(x);
N = size(x,2);
half_bin = size(X,3);

P.Pxii = ones(N,half_bin);
P.Pxij = ones(N,N,half_bin);
P.Fvv = zeros(N,N,half_bin);

G = zeros(1,half_bin);
SNR = zeros(size(X,1),half_bin);

% output spectral
Y = squeeze(X(:,1,:));

alpha = 0.68;
alpha_MSC = 0;
for frameIndex = 1:size(X,1)
    X_t = squeeze(X(frameIndex,:,:));
    P = update_MSC(X_t,P,alpha,alpha_MSC);
    
    Gmin = 0.1;
    method = 4;
    for k = 1:half_bin
        [G(k),SNR(frameIndex,k)] = getweights(P.Fvv,k,d,Gmin,method);
    end
    Y(frameIndex,:) = Y(frameIndex,:).*G;
    
end
y = istft(Y);
Fvv = P.Fvv;

end

