%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% endfire
% refer to "A Dual-Microphone Speech Enhancement Algorithm
% Based on the Coherence Function"
%
% broadside
% refer to "A coherence-based noise reduction algorithm for binaural
% hearing aids"
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clear all;
addpath(genpath('../../../../lib'));
c = 340; % speed of sound

%%
%% load recorded office noise audio

fs = 16000;

angle = [0,0]/180*pi;
% array spacing
d = 0.064;
r = d/2; 

%% generate diffuse noise
[ noise,Pos] = noise_gen_URA(r);

noise = noise'/500;
N = size(noise,2);        %Channels
%% 

pathname = '../../sound/';

%use a clean speech audio as desired signal
[speech ,fs] = audioread([pathname,'S_01_01.wav']);
[interference ,fs] = audioread([pathname,'S_72_09.wav']);

% pathname = '../sound/speech/';
               
%% use RIR to generate reverberation output

angle = [0 0]/180*pi;          % source direction [0,180]
[x1,y1,z1]=sph2cart(angle(1),angle(2),1);    % source position 1
source_pos = [x1,y1,z1];              % Source position [x y z] (m)

angle2 = [90 0]/180*pi;          % source direction [90,180]
[x2,y2,z2]=sph2cart(angle2(1),angle(2),1);    % source position 1
interf_pos2 = [x2,y2,z2];              % Source position [x y z] (m)

%source_pos = [2-d*2 1.5+0.7 2];              % Source position [x y z] (m)

scale = 10;

h = RIR_generator_URA( source_pos,0.2,r);
h = h*scale;
% h = h(:,1:1024);
s = zeros(length(speech)+size(h,2)-1,N);
for i=1:N
    s(:,i) = conv(speech,h(i,:));
end

h = RIR_generator_URA( interf_pos2,0.2,r);
h = h*scale;
% h = h(:,1:1024);
interf = zeros(length(interference)+size(h,2)-1,N);
for i=1:N
    interf(:,i) = conv(interference,h(i,:));
end

% signal+interference+diffuse noise
len_min = min(min(size(s,1),size(interf,1)),size(noise,1));
x = s(1:len_min,:)+interf(1:len_min,:)+noise(1:len_min,:);
x = s(1:len_min,:)+interf(1:len_min,:);
x = interf(1:len_min,:);%+noise(1:len_min,:);
% x = s(1:len_min,:)+noise(1:len_min,:);
% x = s(1:len_min,:);
% x = noise(1:len_min,:);

% x = s(:,1:2);
% x = noise(:,1:2);
x = x(:,[1,3]);
% x = [interference(3:end),interference(1:end-2)];
% d = 0.0213*2;
% x = audioread(['ETAudioDump_245.wav']);
% x = [x(:,2),x(:,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = size(x,2);
%% 
frameLength = 256;
overlap = 128;
inc = frameLength - overlap;
N_FFT = 256;
half_bin = N_FFT/2+1;
c = 340;

tau = fs*d/c;

x1 = x;
[ out,Fvv2] = process(x1,d);
% visual( x(:,1),out );

% x = interf(:,[1,3]);

% Fvv2 = ones(N_FFT/2+1,M,M);     % postfilter MSC
% P_len = N_FFT/2+1;
% Pxii_pre = ones(N,P_len);
% Pxij_pre = ones((N*N-N)/2,P_len);
% 
% x1 = x;
%% estimate noise coherence function
% noise = x1(1:5000,:);
% noise = x1(15000:30000,:);
% noise = x(36000:end,:);
noise = x;
% for i = 1:size(noise,2)
%     for j = 1:size(noise,2)       
%         [sc,F]=mycohere(noise(:,i),noise(:,j),N_FFT,fs,hanning(N_FFT),0.75*N_FFT);
%         Fvv2(:,i,j) = sc;%real(sc);
% 
%     end
% end


k = 0:1:half_bin;
% k = k(2:end);
omega = 2*pi*k/N_FFT;
figure,subplot(2,1,1),
plot(real(squeeze(Fvv2(1,2,:))))
% plot(real(squeeze(Fvv2(:,1,2))))
hold on,plot(cos(omega*tau*cos(angle2(1)))),title('real part'),legend('estimate','theoretic');
subplot(2,1,2),
plot(imag(squeeze(Fvv2(1,2,:))))
hold on,plot(sin(omega*tau*cos(angle2(1)))),title('image part'),legend('estimate','theoretic');

% Fvv = GenNoiseMSC(M,N_FFT,fs,d);  % superdirective noise CSD
% figure,plot(real(squeeze(Fvv2(1,2,:))))
% hold on,plot(Fvv(:,1,2))



