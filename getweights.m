function [ G ] = getweights( Fvv,k,d)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明

alpha_low = 16;
alpha_hi = 2;
beta_low = -0.1;
beta_hi = -0.3;
mu = 0.05;
gamma = 0.6;

SNR_est  = (1-(real(Fvv(1,2,k)))^2-imag(Fvv(1,2,k))^2)/...
           ((1-real(Fvv(1,2,k)))^2+imag(Fvv(1,2,k))^2);
if SNR_est <0.01
    L = 1;
elseif SNR_est>100
    L = 512;
else
    L = 2^(SNR_est/5+5);
end

G = 1-abs(real(Fvv(1,2,k)))^L;

if(k<=16)
    G1 = 1-abs(real(Fvv(1,2,k)))^alpha_low;
    Q = beta_low;
else
    G1 = 1-abs(real(Fvv(1,2,k)))^alpha_hi;
    Q = beta_hi;
end
if(imag(Fvv(1,2,k))<Q)
    G2 = mu;
else
    G2 = 1;
end
G = G1*G2;  % endfire


% G = sqrt(SNR_est/(1+SNR_est));
% 
%                 % broadside
%                 G = sqrt((1-(real(Fvv(1,2,k)))^2-imag(Fvv(1,2,k))^2)/...
%                 (2*(1-real(Fvv(1,2,k)))));

% GAIN_FLOOR = 0.01;
% if(G<GAIN_FLOOR)
%     G = GAIN_FLOOR;%*sign(G(k));
% end
% if(G>=1)
%     G = 1;
% end


end

