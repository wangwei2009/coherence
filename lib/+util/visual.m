function [ output_args ] = visual( x,y )
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明
%% visualization
X = stft(x);
Y = stft(y);
X = abs(squeeze(X)');
Y = abs(squeeze(Y)');
figure,
subplot(211)
imagesc(10*log10(X))
set(gca,'YDir','normal')
caxis([-50 10])
colorbar
title('input signal')
xlabel('frame')
ylabel('frequeny')
subplot(212)
imagesc(10*log10(Y))
set(gca,'YDir','normal')
caxis([-50 10])
colorbar
title('output signal')
xlabel('frame')
ylabel('frequeny')

end

