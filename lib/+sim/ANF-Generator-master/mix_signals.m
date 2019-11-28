function x = mix_signals(n,DC,method)

% Mix M mutually indepedent signals such that the mixed signals 
% exhibit a specific spatial coherence.
%
% Input parameters:
%       n      : M signals in the STFT domain [L x M]
%       DC     : Desired coherence [M x M x K/2+1]
%       method : 'cholesky' or 'eigen'
%
% Output parameters:
%       x      : M generated signals [L x M]
%
% Author       : E.A.P. Habets
% Date         : 29-06-2017
%
% Reference    : E.A.P. Habets, I. Cohen and S. Gannot, 'Generating 
%                nonstationary multisensor signals under a spatial 
%                coherence constraint', Journal of the Acoustical Society
%                of America, Vol. 124, Issue 5, pp. 2911-2917, Nov. 2008.

% Copyright (C) 2009-2017 E.A.P. Habets
%  
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

narginchk(2,3);

if nargin < 3
    method = 'cholesky';
end

M = size(n,2); % Number of sensors
L = size(n,1); % Length input signal
K = (size(DC,3)-1)*2;

% Short-time Fourier transform
for m = 1 : M
    N(m,:,:) = stft(n(:,m), K, K/4, 1).';
end

% Initialization
C = zeros(size(DC)); % STFT mixing matrix
X = zeros(size(N));  % STFT output matrix
X(:,:,1) = X(1,1,1);

% Generate output in the STFT domain for each frequency bin k
for k = 2:K/2+1
    switch lower(method)
        case 'cholesky'
            C(:,:,k) = chol(DC(:,:,k));
            
        case 'eigen'
            [V,D] = eig(DC(:,:,k));
            C(:,:,k) = sqrt(D) * V';
            
        otherwise
            error('Unknown method specified.');
    end
    
    X(:,:,k) = C(:,:,k)' * N(:,:,k);
end

% Inverse STFT
for m = 1 : M
    x(:,m) = real(istft(squeeze(X(m,:,:)).', K, K/4, 1));
end
x = x(1:L,:);