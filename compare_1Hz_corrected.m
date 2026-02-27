%% 修正版：computeBandEnergy vs pwelch 对比（正确积分PSD）
clear; clc; close all;

fprintf('========== 修正版对比测试（正确积分PSD） ==========\n\n');

%% 1. 生成测试信号
Fs = 1000;
T = 10;
t = (0:1/Fs:T-1/Fs)';
N = length(t);

A = 2;
f0 = 100;
signal = A * sin(2*pi*f0*t) + 0.01*randn(N,1);

fprintf('测试信号: A=%.2f Pa, f=%d Hz, N=%d\n\n', A, f0, N);

%% 2. computeBandEnergy（1Hz带宽）
test_freqs = 98:105;
numBands = length(test_freqs);
FreqBands = [test_freqs', test_freqs'+1];

bandPower = computeBandEnergy(signal, Fs, FreqBands);

fprintf('【方法1】computeBandEnergy (1Hz带宽)\n');
for i = 1:numBands
    fprintf('%3d-%3d Hz: %.6f Pa²\n', test_freqs(i), test_freqs(i)+1, bandPower(i));
end
fprintf('\n');

%% 3. pwelch + 正确积分
window = hanning(4096);
noverlap = 2048;
nfft = 16384;
[Pxx, f_psd] = pwelch(signal, window, noverlap, nfft, Fs);

df = f_psd(2) - f_psd(1);  % 频率分辨率
fprintf('【方法2】pwelch + 积分\n');
fprintf('频率分辨率: %.4f Hz\n\n', df);

% 对每个1Hz频带，积分PSD
Power_from_pwelch = zeros(numBands, 1);
for i = 1:numBands
    fLow = test_freqs(i);
    fHigh = test_freqs(i) + 1;

    % 找到该频带内的所有频率点
    idx_band = (f_psd >= fLow) & (f_psd < fHigh);

    % 对PSD积分（求和 × 频率分辨率）
    Power_from_pwelch(i) = sum(Pxx(idx_band)) * df;

    fprintf('%3d-%3d Hz: %.6f Pa² (积分%d个点)\n', ...
            fLow, fHigh, Power_from_pwelch(i), sum(idx_band));
end
fprintf('\n');

%% 4. 对比
fprintf('【对比结果】\n');
fprintf('频带     | computeBandEnergy | pwelch积分 | 差异(dB)\n');
fprintf('---------|-------------------|------------|----------\n');

for i = 1:numBands
    diff_dB = 10*log10(bandPower(i) / Power_from_pwelch(i));
    fprintf('%3d-%3d Hz | %.6f Pa²      | %.6f Pa² | %+6.2f\n', ...
            test_freqs(i), test_freqs(i)+1, bandPower(i), Power_from_pwelch(i), diff_dB);
end

% 重点关注99-100Hz和100-101Hz
idx_99 = find(test_freqs == 99);
idx_100 = find(test_freqs == 100);

fprintf('\n关键频带:\n');
fprintf('  99-100 Hz 差异: %.2f dB\n', 10*log10(bandPower(idx_99)/Power_from_pwelch(idx_99)));
fprintf(' 100-101 Hz 差异: %.2f dB\n', 10*log10(bandPower(idx_100)/Power_from_pwelch(idx_100)));

mean_diff = mean(abs(10*log10(bandPower ./ Power_from_pwelch)));
fprintf('\n平均绝对差异: %.2f dB\n', mean_diff);

%% 5. 绘图
figure('Position', [100, 100, 1200, 800]);

subplot(2,1,1);
bar_data = [bandPower, Power_from_pwelch];
bar(test_freqs + 0.5, bar_data);
legend('computeBandEnergy', 'pwelch积分', 'Location', 'best');
xlabel('频带中心频率 (Hz)');
ylabel('功率 (Pa²)');
title('1Hz带宽功率对比');
grid on;

subplot(2,1,2);
diff_dB_all = 10*log10(bandPower ./ Power_from_pwelch);
bar(test_freqs + 0.5, diff_dB_all);
yline(0, 'r--', 'LineWidth', 2);
xlabel('频带中心频率 (Hz)');
ylabel('差异 (dB)');
title('差异 (computeBandEnergy - pwelch)');
grid on;

%% 6. 结论
fprintf('\n========== 结论 ==========\n');
if mean_diff < 1
    fprintf('✓✓ 当频带宽度为1Hz时，两种方法结果高度一致！\n');
    fprintf('   平均差异 < 1 dB\n');
elseif mean_diff < 3
    fprintf('✓ 当频带宽度为1Hz时，两种方法结果基本一致\n');
    fprintf('   平均差异 < 3 dB\n');
else
    fprintf('○ 存在一定差异（平均 %.2f dB）\n', mean_diff);
    fprintf('   主要由以下原因导致:\n');
    fprintf('   1. computeBandEnergy使用矩形窗（频率泄漏）\n');
    fprintf('   2. pwelch使用Hanning窗（减少泄漏）\n');
    fprintf('   3. pwelch使用分段平均（降低方差）\n');
end
fprintf('==========================\n');
