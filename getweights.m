function [G, SNR] = getweights(Fvv, k, d, Gmin, method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% get dual-mic coherence-based weight
% example Usage:
%   Y = getweights(Fvv,16,0.032,0.1)
%
% Inputs:
%   Fvv           coherence function at k frequency bin
%   k             frequency bin
%   d             dual-min distance
%   Gmin          gain floor
%
%
% Outputs:
%   G            filter gain
%   SNR          estimate SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('Gmin', 'var')
    Gmin = 0.1;
end

if ~exist('method', 'var')
    method = 2;
end

fs = 16000;
N_FFT = (size(Fvv,3)-1)*2;
c = 340;
dij = d;

Fy_real = real(Fvv(1, 2, k));
Fy_imag = imag(Fvv(1, 2, k));

switch method
    case 1
        % endfire
        % refer to 
        % "A Dual-Microphone Speech Enhancement Algorithm
        %  Based on the Coherence Function"
        alpha_low = 16;
        alpha_hi = 2;
        beta_low = -0.1;
        beta_hi = -0.3;
        mu = 0.05;
        gamma = 0.6;

        SNR_est = (1 - (real(Fvv(1, 2, k)))^2 - imag(Fvv(1, 2, k))^2) / ...
            ((1 - real(Fvv(1, 2, k)))^2 + imag(Fvv(1, 2, k))^2);

        if SNR_est < 0.01
            L = 1;
        elseif SNR_est > 100
            L = 512;
        else
            L = 2^(SNR_est / 5 + 5);
        end

        G = 1 - abs(real(Fvv(1, 2, k)))^L;

        if (k <= 16)
            G1 = 1 - abs(real(Fvv(1, 2, k)))^alpha_low;
            Q = beta_low;
        else
            G1 = 1 - abs(real(Fvv(1, 2, k)))^alpha_hi;
            Q = beta_hi;
        end

        if (imag(Fvv(1, 2, k)) < Q)
            G2 = mu;
        else
            G2 = 1;
        end

        G = G1 * G2; % endfire
    case 2
        % refer to 
        % [1] "A Dual-Microphone Algorithm That Can
        %      Cope With Competing-Talker Scenarios"
        % endfire ,coherent

        %theta = 0 * pi / 180;             % 90,interference broadside
        ata = 0 * pi / 180;                % 0,target endfire
        omega = 2 * pi * (k - 1) / N_FFT;
        tao = fs * d / c;
        omega_ = omega * tao;
        beta = omega_ * cos(ata);
        %alpha = omega_ * cos(theta);

        A = Fy_imag - sin(omega_);
        B = cos(omega_) - Fy_real;
        C = Fy_real * sin(omega_) - Fy_imag * cos(omega_);     % eq.13

        T = 1 - Fy_real * cos(omega_) - Fy_imag * sin(omega_); % eq.18

        sin_alpha = (-1 * B * C + A * T) / ...
                    (A^2 + B^2);                               % eq.17/18
        SNR = (sin_alpha - Fy_imag) / ...
              (Fy_imag - sin(beta));                           % eq.11 in [1]

        % square-root wiener filter
        G = sqrt(SNR / ...
            (SNR + 1));                           

    case 5
        % refer to 
        % [1] "A coherence-based noise reduction algorithm for binaural
        %      hearing aids"
        % broadside ,coherent
        % [2] "A Joint Speech Enhancement Algorithm Based on the
        %      Tri-microphone"

        sin_alpha = (2*(1-Fy_real)*Fy_imag) / ...
                     ((1-Fy_real)^2+Fy_imag^2);          % eq.19 in [1]
        G = (1-Fy_real^2-Fy_imag^2)/...
            (2*(1-Fy_real));                             % eq.20 in [1]
        SNR = (1-Fy_real^2-Fy_imag^2)/...
              ((1-Fy_real)^2+Fy_imag^2);                 % eq.23 in [1]
                                                         %    is
                                                         %   equal
        SNR = (sin_alpha - Fy_imag) / ...                %    to
              (Fy_imag);                                 % eq.7 in [2]
        G = sqrt(SNR / ...
            (SNR + 1));
    case 3
        %
        % refer to 
        % "Hybrid Coherence Model for Noise Reduction in Reverberant
        %  Environments"
        %
        k_optimal = 1;


        abs_Fvv2 = sqrt(Fy_real^2 + Fy_imag^2);
        Fn = sin(2 * pi * k * fs * dij * k_optimal / c / N_FFT) ./ (2 * pi * k * fs * dij * k_optimal / c / N_FFT);
        %                 Fn = sinc(2*pi*k*fs*d/(N_FFT*c));

        DDR = (abs(Fn)^2 - abs_Fvv2^2) / ...
            (abs_Fvv2^2 - 1);
        K = DDR / (DDR + 1);
        %theta = 0 * pi / 180; % 90,interference broadside
        ata = 0 * pi / 180; % 0,target endfire
        %                 omega = 2*pi*k/N_FFT;
        omega = 2 * (k - 1) / N_FFT;
        tao = fs * d / c;
        omega_ = omega * tao;
        beta = pi * omega_ * cos(ata);
        %alpha = omega_ * cos(theta);
        constant = 2 * pi * k * fs * d / ((N_FFT * c));

        % if we set K = 1,then this method is same as method2
        A = K * (Fy_imag - sin(beta));
        B = (1 - K) * sinc(omega_) + K * cos(beta) - Fy_real;
        C = Fy_real * sin(beta) - Fy_imag * K * cos(beta) - (1 - K) * sinc(omega_) * sin(beta);
        T = K - Fy_real * cos(beta) - K * Fy_imag * sin(beta) + (1 - K) * sinc(omega_) * cos(beta);

        sin_alpha = (-1 * B * C + A * T) / ...
            (A^2 + B^2); % eq.14
        SNR = (sin_alpha - Fy_imag) / ...
            (Fy_imag - sin(beta)); % eq.10

        if (SNR <- 0.99 || isnan(SNR))
            SNR = -0.99;
        end
        G = SNR / ...
            SNR + 1;
        G = G^2;

%             G = sqrt(SNR / ...
%                 (SNR + 1));
    case 4
        % refer to "Coherence-based dual-channel noise reduction algorithm in a complex noisy
        % environment" method_1
        % endfire ,coherent+diffuse
        k_optimal = 1;

        Fy_real = real(Fvv(1, 2, k));
        Fy_imag = imag(Fvv(1, 2, k));
        Fn = sin(2 * pi * k * fs * dij * k_optimal / c / N_FFT) ./ (2 * pi * k * fs * dij * k_optimal / c / N_FFT);
        %  Fn = sinc(2*pi*k*fs*d/(N_FFT*c));


        abs_Fvv2 = sqrt(Fy_real^2 + Fy_imag^2);


        DDR = (abs(Fn)^2 - abs_Fvv2^2) / ...
            (abs_Fvv2^2 - 1+1e-6);
%       DDR = max(0,DDR);
        K = DDR / (DDR + 1);
        %                 K = 1;
        theta_s = 90 * pi / 180; % 90,target,endfire
        %theta_i = 0 * pi / 180; % 0,interference ,broadside
        constant = 2 * pi * (k - 1) * fs * d / ((N_FFT * c));
        sin_alpha = sin(constant * sin(theta_s));
        cos_alpha = cos(constant * sin(theta_s));
        A = sin_alpha * K - Fy_imag;
        B = cos_alpha * K - Fy_real + Fn * (1 - K);
        C = (Fy_real - Fn * (1 - K)) * sin_alpha - Fy_imag * cos_alpha;
        T = K - cos_alpha * (Fy_real - Fn * (1 - K)) - Fy_imag * sin_alpha;
        sin_beta = (-1 * B * C - A * T) / ...
            (A^2 + B^2);
        G = ((Fy_imag - sin_beta * K) / ...
            (sin_alpha - sin_beta));
        SNR = G / (1 - G);
%             G = G*K;
    case 6
        % refer to "Coherence-based dual-channel noise reduction algorithm in a complex noisy
        % environment" method_3
        % endfire ,coherent+diffuse
        k_optimal = 1;
        SPECTRAL_FLOOR = 0.4;
        Fvv_UPPER = 0.98;
        Fy_real = real(Fvv(1, 2, k));
        Fy_imag = imag(Fvv(1, 2, k));
        Fn = sin(2 * pi * k * fs * dij * k_optimal / c / N_FFT) ./ (2 * pi * k * fs * dij * k_optimal / c / N_FFT);
        %                 Fn = sinc(2*pi*k*fs*d/(N_FFT*c));

        if (Fy_real > Fvv_UPPER)
            Fy_real = Fvv_UPPER;
        end

        abs_Fvv2 = sqrt(Fy_real^2 + Fy_imag^2);

        if (abs_Fvv2 > Fvv_UPPER)
            abs_Fvv2 = Fvv_UPPER;
        end


        DDR = (abs(Fn)^2 - abs_Fvv2^2) / ...
            (abs_Fvv2^2 - 1);
        K = DDR / (DDR + 1);
        %                 K = 1;
        theta_s = 90 * pi / 180; % 90,target,endfire
        %theta_i = 0 * pi / 180; % interference ,broadside
        constant = 2 * pi * (k - 1) * fs * d / ((N_FFT * c));
        sin_alpha = sin(constant * sin(theta_s));
        cos_alpha = cos(constant * sin(theta_s));

        A = sin_alpha * K - Fy_imag;
        B = cos_alpha * K - Fy_real + Fn * (1 - K);
        C = (Fy_real - Fn * (1 - K)) * sin_alpha - Fy_imag * cos_alpha;
        T = K - cos_alpha * (Fy_real - Fn * (1 - K)) - Fy_imag * sin_alpha;
        sin_beta = (-1 * B * C - A * T) / ...
                   (A^2 + B^2);                % eq.21
        cos_beta = (A*C-B*T)/(A^2+B^2);        % eq.22

        A_ = cos_alpha - cos_beta;
        B_ = cos_beta + Fn*(1-K);
        C_ = sin_alpha - sin_beta;
        D_ = sin_beta*K;              % eq.16 

        if(abs(Fy_imag-sin_alpha)<abs(Fy_imag-sin_beta))
            gamma = 1;
        else
            gamma = -1;
        end                           % eq.18

        T_ = abs_Fvv2^2*(A_^2+C_^2)-(A_*D_-B_*C_)^2;  

        G = (-1*(A_*B_+C_*D_)+gamma*sqrt(T_))/(A_^2+C_^2);   % eq.17

        SNR = G / (1 - G);
    case 7
        % refer to "Roubust Recognition of Reverberant and noisy speech using 
        % coherence-based processing" 
        % braodside ,coherent+diffuse
        k_optimal = 1;

        Fn = sin(2 * pi * k * fs * dij * k_optimal / c / N_FFT) ./ (2 * pi * k * fs * dij * k_optimal / c / N_FFT);
        % Fn = sinc(2*pi*k*fs*d/(N_FFT*c));

        DDR = (Fn - Fy_real) / ...
            (Fy_real - 1+1e-6);
        DDR = max(0,DDR);
        K = DDR / (DDR + 1);

        SNR = K;
        G = K;
end

G = max(G, Gmin);
G = min(G, 1);

end
