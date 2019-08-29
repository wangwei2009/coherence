%SPECTRAL_SUBTRACTION Compute spectral subtraction weights.
%
% weights = spectral_subtraction(SNR,alpha,beta,mu,Gmin)
%
% alpha = 1; beta = 1;   % power subtraction
% alpha = 2; beta = 0.5; % magnitude subtraction
% alpha = 2; beta = 1;   % Wiener filter
% mu: noise overestimation
% Gmin: gain floor
%
% Andreas Schwarz (schwarz@lnt.de)
% Multimedia Communications and Signal Processing
% Friedrich-Alexander-Universitaet Erlangen-Nuernberg (FAU)
% Cauerstr. 7, 91058 Erlangen, Germany
function weights = spectral_subtraction(SNR,alpha,beta,mu,Gmin)

if (nargin == 1)
    % default: magnitude subtraction
    alpha = 2;
    beta = 0.5;
    mu = 1;
end

if ~exist('Gmin','var')
    Gmin = 0.1;
end

SNR = max(SNR,0);
weights = max(1 - (mu./(SNR + 1)).^beta, 0).^alpha;
weights = max(weights,0);
weights = max(sqrt(weights),Gmin);

end