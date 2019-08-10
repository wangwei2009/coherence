function [ Fvv ] = GenNoiseMSC(M,N_FFT,fs,r,varargin)
%GenNoiseMSC compute diffuse noise field Mignitude-SquareD Coherence
%function
%
% Usage: (1) [ Fvv ] = GenNoiseMSC(M,N,fs,r,arrayType);  
%
% Inputs:
%   M              input channels
%   N_FFT          fft points
%   fs             sample frequency in Hz
%   r              array aperture
%   varargin
%       arrayType  array type,linear or circular,default linear
%
% Outputs:
%   Fvv            [half_bin,M,M] coherence function
if(nargin>4)
    arrayType = varargin{1};
else
    arrayType = 'linear';
end
c = 340;
f = 0:fs/256:fs/2;
Fvv = zeros(N_FFT/2+1,M,M);
k_optimal = 1;
for i = 1:M
    for j = 1:M   
        if i == j
            Fvv(:,i,j) = ones(N_FFT/2+1,1);
        else
            switch(lower(arrayType))
                case 'linear'
                    dij = abs(i-j)*r;
                case 'circular'
                    mic_rad = abs(i-j)*(360/M)*pi/180;         % radian between two sensor
                    dij = np.sqrt(r^2+r^2-2*r*r*cos(mic_rad)); % distant between two sensor
            end
            Fvv(:,i,j) = sin(2*pi*f*dij*k_optimal/c)./(2*pi*f*dij*k_optimal/c);Fvv(1,i,j) = 0.998;%T(2) = 0.996;
        end
%         index = find(Fvv(:,i,j)>0.9);
%         if(size(index,1)>0)
%             Fvv(index,i,j)=0.9;
%         end
    end
end

end

