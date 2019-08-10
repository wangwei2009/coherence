function [ G ] = getweights3(Fvv,k,d)
%
%
% refer to "Hybrid Coherence Model for Noise Reduction in Reverberant
% Environments"
%
%
fs = 16000;
dij = d;
k_optimal = 1;
c = 340;
N_FFT = 256;
                
                Fvv_UPPER = 0.98;
                Fy_real = real(Fvv(1,2,k));
                Fy_imag = imag(Fvv(1,2,k));
                if(Fy_real>Fvv_UPPER)
                    Fy_real = Fvv_UPPER;
                end
                abs_Fvv2 = sqrt(Fy_real^2+Fy_imag^2);
                Fn = sin(2*pi*k*fs*dij*k_optimal/c/N_FFT)./(2*pi*k*fs*dij*k_optimal/c/N_FFT);
%                 Fn = sinc(2*pi*k*fs*d/(N_FFT*c));
                if(Fn>Fvv_UPPER)
                    Fn = Fvv_UPPER;
                end
                
                DDR = (abs(Fn)^2-abs_Fvv2^2)/...
                      (abs_Fvv2^2-1);
                K = DDR/(DDR+1);
                theta = 90*pi/180; % interference broadside
                ata = 0*pi/180;      % target endfire  
%                 omega = 2*pi*k/N_FFT;
                omega = 2*(k-1)/N_FFT;
                tao = fs*d/c;
                omega_ = omega*tao;
                beta = pi*omega_*cos(ata);
                alpha = omega_*cos(theta);
                constant = 2*pi*k*fs*d/((N_FFT*c));
                
%                 K = 1;
                A = K*(Fy_imag-sin(beta));
                B = (1-K)*sinc(omega_)+K*cos(beta)-Fy_real;
                C = Fy_real*sin(beta)-Fy_imag*K*cos(beta)-(1-K)*sinc(omega_)*sin(beta);
                T = K-Fy_real*cos(beta)-K*Fy_imag*sin(beta)+(1-K)*sinc(omega_)*cos(beta);

                sin_alpha = (-1*B*C+A*T)/...
                           (A^2+B^2);               % eq.14
                SNR = (sin_alpha-Fy_imag)/...
                      (Fy_imag-sin(beta));          % eq.10
                if(SNR<-0.99||isnan(SNR))
                    SNR = -0.99;
                end
                G = spectral_subtraction(SNR,2,0.5,10);
                G = sqrt(SNR/...
                           (SNR+1));

                GAIN_FLOOR = 0.1;
                if(G<GAIN_FLOOR)
                    G = GAIN_FLOOR;%*sign(G(k));
                end
                if(G>=1)
                    G = 1;
                end
                if(isnan(G))
                    G = GAIN_FLOOR;
                end
                    
            end
                

