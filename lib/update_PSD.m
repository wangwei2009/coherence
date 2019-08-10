function [ Pxii ] = update_PSD( X,alpha,Pxii)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate auto-spectral power density
% 
%
% example Usage: 
%   [ Pxii ] = update_PSD( X,alpha,Pxii)
%
% Inputs:
%   X              Multichannel input data,stored in column-wise
%   alpha          average factor
%   Pxii           last Pxii,size of [N(N-1)/2,half_bin]
%       
%
% Outputs:
%   Pxii            recursively averaged PSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(X,1); % channel
for i = 1:N
    Xi = X(i,:);
    % update auto-spectral power density
    Pxii(i,:) = alpha*Pxii(i,:)+(1-alpha)*Xi.*conj(Xi);
end

end

