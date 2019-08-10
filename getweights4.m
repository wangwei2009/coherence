function [ G,SNR ] = getweights4( Fvv,k,d )
% refer to "A Dual-Microphone Algorithm That Can
%           Cope With Competing-Talker Scenarios"
% endfire ,coherent
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
         
                DDR = (abs(Fn)^2-abs(Fvv(1,2,k))^2)/...
                      (abs(Fvv(1,2,k))^2-1);
                K = DDR/(DDR+1);
                theta = 90*pi/180; % interference broadside
                ata = 0*pi/180;      % target endfire  
                omega = 2*pi*(k-1)/N_FFT;
                tao = fs*d/c;
                omega_ = omega*tao;
                beta = omega_*cos(ata);
                alpha = omega_*cos(theta);
                constant = 2*pi*k*fs*d/((N_FFT*c));

%                 K = 1;
                A = Fy_imag-sin(omega_);
                B = cos(omega_)-Fy_real;
                C = Fy_real*sin(omega_)-Fy_imag*cos(omega_);
                T = 1-Fy_real*cos(omega_)-Fy_imag*sin(omega_);

                sin_alpha = (-1*B*C+A*T)/...
                           (A^2+B^2);               % eq.14
                SNR = (sin_alpha-Fy_imag)/...
                      (Fy_imag-sin(beta));          % eq.10
                G = sqrt(SNR/...
                           (SNR+1));
                       
                GAIN_FLOOR = 0.01;
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

