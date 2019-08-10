function [ h ] = RIR_generator_URA( s,beta,r)
%UNTITLED7 此处显示有关此函数的摘要
%   此处显示详细说明
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
d = 0.05;
orig = [2 1.5 3];
s = s+orig;
% r = 0.032;
[x1,y1,z1]=sph2cart(0,0,r);    % Sensor position 1
[x2,y2,z2]=sph2cart(90/180*pi,0,r);    % Sensor position 2
[x3,y3,z3]=sph2cart(180/180*pi,0,r);    % Sensor position 1
[x4,y4,z4]=sph2cart(270/180*pi,0,r);    % Sensor position 2
P = [x1 x2 x3 x4;y1 y2 y3 y4;z1 z2 z3 z4]; % Construct position matrix

r = [2   1.5 2 ;
    2-d 1.5 2;
    2-d*2 1.5 2 ;
    2-d*3 1.5 2;
    2-d*4 1.5 2;];
r = ones(4,3).*orig+P';% Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
%s = [1.94 3.5 2];              % Source position [x y z] (m)
L = [5 4 6];                % Room dimensions [x y z] (m)
%beta = 0.2;                 % Reverberation time (s)
n = beta*fs;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 0;              % Enable high-pass filter

h = rir_generator(c, fs, r, s, L, beta,n, mtype, order, dim, orientation, hp_filter);

end

