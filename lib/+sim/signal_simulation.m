function [ sig ] = signal_simulation( r,slice )
% signal simulation function
%
% Usage: (1) [ sig ] = signal_simulation(0.064,[1,3]); 
%
% Inputs:
%   r              radius of circular array
%   slice          choose which channel to use
%   varargin
%       arrayType  array type,linear or circular,default linear
%
% Outputs:
%   sig            4-channel siganl [samples,channel]
%
% Created by Wang wei
%% generate diffuse noise
% [ noise,Pos] = sim.noise_gen_URA(r);
addpath(genpath('ANF-Generator-master'));
noise = sim.gen_babble_speech(4,r,5);
% noise = noise'/500;
%% signal simulation 
pathname = 'wav/';

% use a clean speech audio as desired signal
[speech ,fs] = audioread([pathname,'S_01_01.wav']);
[interference] = audioread([pathname,'S_72_09.wav']);

s      = sim.generate_signal(speech,[0,0],r,1,0.2);
interf = sim.generate_signal(interference,[90,0],r,1,0.2);

% signal+interference+diffuse noise
len_min = min(min(size(s,1),size(interf,1)),size(noise,1));
% slice = [1,3];
clean_s = s(1:len_min,slice);
clean_i = interf(1:len_min,slice);
clean_n = noise(1:len_min,slice)/2;
switch 1
    case 1
        x = clean_s+clean_i+clean_n;
    case 2
        x = clean_s+clean_i;
        clean_n = zeros(size(clean_n));
    case 3
        x = clean_s;
        clean_i = zeros(size(clean_i));
        clean_n = zeros(size(clean_n));
    case 4
        % ideal signal
        len_min_clean = min(length(speech),length(interference));
        speech = speech(1:len_min_clean);
        interference = interference(1:len_min_clean);
        clean_s = [speech(3:end),speech(1:end-2)];
        clean_i = [interference(1:size(clean_s,1)),interference(1:size(clean_s,1))];
        x = clean_s + clean_i;
        d = 0.0213*2;
    otherwise
        disp('other value')
end

SNR = snr(clean_s(:,1),clean_i(:,1))

sig.speech = speech;
sig.clean_s = clean_s;
sig.clean_i = clean_i;
sig.clean_n = clean_n;
sig.x = x;
rmpath(genpath('./ANF-Generator-master'));
end

