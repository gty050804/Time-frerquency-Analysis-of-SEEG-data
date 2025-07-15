% 生成一个非平稳信号
fs = 1000;          % 采样频率 (Hz)
t = 0:1/fs:2;       % 时间向量
x = cos(2*pi*50*t) + cos(2*pi*120*t).*(t > 0.5 & t < 1.5); % 混合信号

% 使用小波变换计算信号的时频表示
[cwt_result, f] = cwt(x, fs);

% 绘制小波变换的时频图
figure;
surface(t, f, abs(cwt_result), 'EdgeColor', 'none');
axis tight;
view(0, 90);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('CWT Magnitude');
colorbar;


cwt_result = abs(cwt_result);