function [ G,SNR] = getweights2(  Fvv,k,d)
% refer to "Coherence-based dual-channel noise reduction algorithm in a complex noisy
% environment"
% endfire ,coherent+diffuse
fs = 16000;
N_FFT = 256;
c = 340;
dij = d;
k_optimal = 1;
SPECTRAL_FLOOR = 0.4;
Fvv_UPPER = 0.98;
                Fy_real = real(Fvv(1,2,k));
                Fy_imag = imag(Fvv(1,2,k));
                Fn = sin(2*pi*k*fs*dij*k_optimal/c/N_FFT)./(2*pi*k*fs*dij*k_optimal/c/N_FFT);
%                 Fn = sinc(2*pi*k*fs*d/(N_FFT*c));
                
                if(Fy_real>Fvv_UPPER)
                    Fy_real = Fvv_UPPER;
                end
                abs_Fvv2 = sqrt(Fy_real^2+Fy_imag^2);
                if(abs_Fvv2>Fvv_UPPER)
                    abs_Fvv2=Fvv_UPPER;
                end
                if(Fn>Fvv_UPPER)
                    Fn = Fvv_UPPER;
                end
                
                
                DDR = (abs(Fn)^2-abs_Fvv2^2)/...
                      (abs_Fvv2^2-1);
                K = DDR/(DDR+1);
%                 K = 1;
                theta_s = 90*pi/180;    % target,endfire
                theta_i = 0*pi/180;     % interference ,broadside
                constant = 2*pi*(k-1)*fs*d/((N_FFT*c));
                sin_alpha = sin(constant*sin(theta_s));
                cos_alpha = cos(constant*sin(theta_s));
                A = sin_alpha*K-Fy_imag;
                B = cos_alpha*K-Fy_real+Fn*(1-K);
                C = (Fy_real-Fn*(1-K))*sin_alpha-Fy_imag*cos_alpha;
                T = K - cos_alpha*(Fy_real-Fn*(1-K))-Fy_imag*sin_alpha;
                sin_beta = (-1*B*C-A*T)/...
                           (A^2+B^2);
                G = ((Fy_imag-sin_beta*K)/...
                       (sin_alpha-sin_beta));
%                (Fy_imag-sin_beta*K)/(sin_alpha-sin_beta)
            
                
%                 % apply CDR estimator (=SNR)
%                 SNR = max(real(SNR),0);
% 
%                 G = spectral_subtraction(SNR,2,0.5,1);
%                 weights = max(weights,0.1);
%                 weights = min(weights,1);
% 
%                 % postfilter input is computed from averaged PSDs of both microphones
%                 Postfilter_input = sqrt((abs((Xi+Xj)/2).^2)) .* exp(1j*angle(Xi));
% 
%                 % apply postfilter
%                 Processed = weights .* Postfilter_input;

                GAIN_FLOOR = 0.1;
                if(G<GAIN_FLOOR)
                    G = GAIN_FLOOR;%*sign(G(k));
                end
                if(G>=1)
                    G = 0.99999;
                end
                if(isnan(G))
                    G = GAIN_FLOOR;
                end
                
                SNR = G/(1-G);
%                 SNR = G^2/(1-G^2);

end

