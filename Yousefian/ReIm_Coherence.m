function [enhanced_ouput]=ReIm_Coherence(x1,x2,fs,d,output_fn)
% x1 , x2 : Input signals at two channels
% fs : Sampling Frequency
% output_path: address to write enhanced file
% Nima Yousefian , Sep 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference: N. Yousefian and P. C. Loizou, "A Dual-Microphone Algorithm 
% for Competing Talker Scenarios", IEEE Transactions on Audio, Speech, 
% and Language Processing, Vol. 21, No. 1, pp. 145-155, Jan. 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Run Speech Enhacement ReIm');
frameLength=floor(15*fs/1000);  
if (frameLength/2 ~= floor(frameLength/2))
      frameLength=frameLength+1;
end
 
frameShift=floor(frameLength* 0.5);
window=hanning(frameLength);
FFT_LEN=2^nextpow2(frameLength);
  
lambda_x=0.68;   %Forgetting factor for smoothing power spectrum 
epsilon=10^-12;
    
lenS=min(length(x1),length(x2));
nFrame=0;
iniFrameSample=1;
endFrameSample=iniFrameSample+frameLength-1;
enhanced_ouput=zeros(lenS,1);
TotFrameNum=floor(lenS/frameShift);
  
dis=d;%0.015;          %Microphones distance (m)
c=340;              %Sound Speed (m/s)
tau=fs * (dis /c);      
Target_DOA= 0 *pi/180;    %DOA -target angle (0 -> Endfire, 90 -> Broadside) 
     
%setting freq. bins between 0 and pi/2
f_bin=0:1/FFT_LEN:0.5;     
f_bin(1)=[];
omega=2 * pi * f_bin;
omega_tau=omega'*tau.*cos(Target_DOA);  
pri_min=10^(-40/10);
  
while endFrameSample<lenS
        
        %A new frame to process
        nFrame=nFrame+1;
      
        %Get short-time magnitude and phase spectrums for each input channel
        Frame1=x1(iniFrameSample:endFrameSample);
        Frame2=x2(iniFrameSample:endFrameSample);
        wFrame1=Frame1 .* window;
        wFrame2=Frame2 .* window;
        X1=fft(wFrame1,FFT_LEN);
        X2=fft(wFrame2,FFT_LEN);
        
        if (nFrame==1)
            PX1X1=abs(X1).^2;
            PX2X2=abs(X2).^2;
            PX1X2=X1.*conj(X2); 
        else
            %auto and cross PSD estimation
            PX1X1=lambda_x.*PX1X1+(1-lambda_x).*abs(X1).^2;
            PX2X2=lambda_x.*PX2X2+(1-lambda_x).*abs(X2).^2;
            PX1X2=lambda_x.*PX1X2+(1-lambda_x).*X1.*conj(X2);            
        end
        
        cohX=PX1X2 ./ sqrt(PX1X1.*PX2X2+epsilon); 
        reCOH=real(cohX(2:FFT_LEN/2+1));
        imCOH=imag(cohX(2:FFT_LEN/2+1));
            
       %Known Parameters - Eq (13)
       A=imCOH-sin(omega_tau);
       B=cos(omega_tau)-reCOH;
       C=reCOH.*sin(omega_tau)- imCOH .* cos(omega_tau);
  
       %Solving unknown parameter -> Alpha
       T=1-reCOH.*cos(omega_tau)-imCOH.*sin(omega_tau);       
      
       sin_alpha=(-B.*C + A .* T)./ (A.^2 + B.^2);        
       SNR_Hat= (sin_alpha - imCOH) ./ (A)  ; 
       SNR_Hat=max(SNR_Hat,pri_min);
         
             %Wiener gain
       W_cur = SNR_Hat(1:FFT_LEN/2) ./ (SNR_Hat(1:FFT_LEN/2) +1);       
       G= sqrt(W_cur);
       H=abs([G ;flipud(G)]);    %Fullband final filter 
    
       %IFFT and OLA
       enhSpeech_Frame_tmp=real(ifft( H .* X1,FFT_LEN));        
       enhSpeech_Frame=enhSpeech_Frame_tmp(1:frameLength);
       enhanced_ouput(iniFrameSample:endFrameSample)=enhSpeech_Frame + enhanced_ouput(iniFrameSample:endFrameSample);      
         
        %Update frame boundaries
        iniFrameSample=iniFrameSample+frameShift;
        endFrameSample=endFrameSample+frameShift;
end
    disp('Processing ... 100% Done');
%     wavwrite(enhanced_ouput,fs,16,output_fn);