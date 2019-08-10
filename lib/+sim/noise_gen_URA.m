function [ z,P ] = noise_gen_URA( r )
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
%% Initialization
fs = 16000;                       % Sample frequency
NFFT = 512;                      % Number of frequency bins (for analysis)
w = 2*pi*fs*(0:NFFT/2)/NFFT; 
c = 340;                         % Speed of sound
L = 2^18;                        % Data length
% r = 0.032;
[x1,y1,z1]=sph2cart(0,0,r);    % Sensor position 1
[x2,y2,z2]=sph2cart(90/180*pi,0,r);    % Sensor position 2
[x3,y3,z3]=sph2cart(180/180*pi,0,r);    % Sensor position 1
[x4,y4,z4]=sph2cart(270/180*pi,0,r);    % Sensor position 2
P = [x1 x2 x3 x4;y1 y2 y3 y4;z1 z2 z3 z4]; % Construct position matrix
M = size(P,2);                           % Number of sensors

% Calculate sensor distances w.r.t. sensor 1
d = zeros(1,M);
for m = 2:M
    d(m) = norm(P(:,m)-P(:,1),2);
end

%% Generate sensor signals
params.c = c;
params.fs = fs;

% 1D example
params.N_phi = 64;
z = sinf_1D(d,L,params); 

% 3D example
% params.N = 256;
% z = sinf_3D(P,L,params); 


end

