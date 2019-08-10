function [ Fvv ] = GenNoiseMSC_shift(M,N,fs,r,DS_weight)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

N_FFT = N;
c = 340;
f = 0:fs/N_FFT:fs/2;
Fvv = zeros(N_FFT/2+1,M,M);
k_optimal = 1;
for i = 1:M
    for j = 1:M   
        if i == j
            Fvv(:,i,j) = ones(N_FFT/2+1,1);
        else
            if(abs(i-j)==1)
                dij = r*sqrt(2);
            elseif(abs(i-j)==2)
                dij = r*2;
            elseif(abs(i-j)==3)
                dij = r*sqrt(2);
            end
            Fvv(:,i,j) = sin(2*pi*f*dij*k_optimal/c)./(2*pi*f*dij*k_optimal/c);Fvv(1,i,j) = 0.998;%T(2) = 0.996;
            Fvv(:,i,j) = Fvv(:,i,j).*DS_weight(:,1).*conj(DS_weight(:,2));
        end
        index = find(Fvv(:,i,j)>0.8);
        if(size(index,1)>0)
            Fvv(index,i,j)=0.8;
        end
    end
end
Fvv = real(Fvv);
Fvv = permute(Fvv,[2,3,1]);
end

