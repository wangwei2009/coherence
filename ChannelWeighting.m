function [aad_X_tilde] = ChannelWeighting(ad_X_Bar,ad_mu,aad_H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gammatone channel weighting
% example Usage:
%   Y = ChannelWeighting(Y,G,H)
%
% Inputs:
%   ad_X_Bar          input frequency dignal
%   ad_mu             weights which to be smoothed
%   aad_H             Gammatone filter
%
% Outputs:
%   aad_X_tilde       smoothed frequency signal
%
% Created By Wang Wei
% refer to 
% "Signal Separation for Robust Speech Recognition Based on
%  Phase Difference Information Obtained in the Frequency Domain"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if ~exist('iNumFilts', 'var')
%     iNumFilts = 40;
% end

iNumFilts = size(aad_H,2);

dEta         = 1e-2;

iFFTSize = (length(ad_X_Bar)-1)*2;

aad_X_tilde = zeros(iFFTSize/2,1);
aad_mu_g    = zeros(iFFTSize / 2,1);
    bPDCW = 1;
if bPDCW == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Calculation of w(i, m)
    %
    for j = 1 : iNumFilts
        aad_w(j) = sum(ad_mu(1 : iFFTSize / 2)' .*   abs (aad_H( : , j)  .*  ad_X_Bar(1 : iFFTSize / 2)').^ 1) ...
                                     / sum(abs(aad_H( :, j) .* ad_X_Bar(1 : iFFTSize / 2)').^ 1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Calculation of /mu_g(k, m) (time-channel mask)
    %
    for j = 1 : iNumFilts,
        aad_mu_g(:) = aad_mu_g(:) + aad_w(j) .^ 1 .* abs(aad_H(:, j)) .^ 1;
    end
    aad_mu_g(:) = ((aad_mu_g(:)) ...
                                   ./ sum(abs(aad_H' .^ 1))') .^ (1 / 1); 

    % Make the mu_g symmetric
%     ad_mu_g_sym = [aad_mu_g(:, iFI); flipud(aad_mu_g( :, iFI))];
%     ad_mu_g_sym = max(ad_mu_g_sym, dEta);
%     aad_X_tilde( :, iFI) = (ad_X_Bar .*  ad_mu_g_sym')';
    aad_mu_g = max(aad_mu_g, dEta);
    aad_X_tilde = ad_X_Bar(1 : iFFTSize / 2) .*  aad_mu_g';
    aad_X_tilde = [aad_X_tilde,0];
end
end

