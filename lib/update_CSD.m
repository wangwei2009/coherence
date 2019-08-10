function [ Pxij ] = update_CSD( X,alpha,Pxij)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate cross-spectral power density
% 
%
% example Usage: 
%   [ Pxij ] = update_CSD( X,alpha,Pxij)
%
% Inputs:
%   X              Multichannel input data,stored in column-wise
%   alpha          average factor
%   Pxij           last Pxij,size of [N(N-1)/2,half_bin]
%       
%
% Outputs:
%   Pxij            recursively averaged CSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(X,1); % channel
for i = 1:N
    for j = 1:N
        Xi = X(i,:);
        Xj = X(j,:);
        if(i==j)
            Fvv(i,j,:) = 1;
        else           
            Pxij_i_j = squeeze(Pxij(i,j,:)).';
            Pxij(i,j,:) = alpha*Pxij_i_j+(1-alpha)*Xi.*conj(Xj);
            
        end
    end
end

