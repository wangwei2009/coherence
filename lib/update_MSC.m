function [ P ] = update_MSC(X,P,alpha,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate coherence function
% 
%
% example Usage: 
%   [ Pxii,Pxij,Fvv ] = update_MSC(X,Pxii,Pxij,Fvv,alpha)
%
% Inputs:
%   X              Multichannel input data,size of [channel,half_bin]
%   alpha          average factor
%   Pxii           previous Pxii,size of [N(N-1)/2,half_bin]
%   Pxij           previous Pxij,size of [N(N-1)/2,half_bin]
%   Fvv            previous Fvv,size of [N,N,,half_bin]
%   varargin
%     alpha_MSC    MSC average factor
%       
%
% Outputs:
%   Pxii,Pxij,Fvv  recursively averaged Pxii,Pxij,Fvv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(X,1); % channel


if nargin <3
    alpha = 0.6;
end
if nargin ==4
    alpha_MSC = varargin{1};
else
    alpha_MSC = 0;
end

Pxii = P.Pxii;
Pxij = P.Pxij;
Fvv = P.Fvv;

% Pxii = update_PSD(X,alpha,Pxii);
% Pxij = update_CSD(X,alpha,Pxij);

for i = 1:N
    Xi = X(i,:);
    % update auto-spectral power density
    Pxii(i,:) = alpha*Pxii(i,:)+(1-alpha)*Xi.*conj(Xi);
end
for i = 1:N
    for j = 1:N
        Xi = X(i,:);
        Xj = X(j,:);
        if(i==j)
            Fvv(i,j,:) = 1;
        else           
            Pxij_i_j = squeeze(Pxij(i,j,:)).';
            Pxij(i,j,:) = alpha*Pxij_i_j+(1-alpha)*Xi.*conj(Xj);
            Pxij_i_j = squeeze(Pxij(i,j,:)).';
            
            % complex coherence function
            Fvv_i_j_curr = Pxij_i_j./(sqrt(Pxii(i,:).*Pxii(j,:)));
            
            Fvv_i_j = squeeze(Fvv(i,j,:)).';
            Fvv(i,j,:) = alpha_MSC*Fvv_i_j+(1-alpha_MSC)*Fvv_i_j_curr;
%             Fvv_i_j_t = Pxij(t,:)./(sqrt(Pxii(i,:).*Pxii(j,:)));
        end
    end
end

P.Pxii = Pxii;
P.Pxij = Pxij;
P.Fvv = Fvv;

