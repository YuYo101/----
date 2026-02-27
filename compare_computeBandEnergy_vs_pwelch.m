%% 对比 computeBandEnergy 与 pwelch 的结果（1Hz带宽）
% 验证：当频带宽度为1Hz时，computeBandEnergy 与 pwelch 的结果是否一致
clear; clc; close all;

fprintf('========== computeBandEnergy vs pwelch 对比测试 ==========\n\n');

%% 1. 生成测试信号
Fs = 1000;  % 1 kHz 采样频率
T = 10;     % 10 秒（较长信号以提高频率分辨率）
t = (0:1/Fs:T-1/Fs)';
N = length(t);

% 单频正弦波：幅度 A = 2 Pa，频率 100 Hz
A = 2;
f0 = 100;
signal = A * sin(2*pi*f0*t);

% 添加少量噪声
signal = signal + 0.01 * randn(N, 1);

fprintf('测试信号:\n');
fprintf('  幅度 A = %.2f Pa\n', A);
fprintf('  频率 f = %d Hz\n', f0);
fprintf('  采样频率 = %d Hz\n', Fs);
fprintf('  信号时长 = %.1f 秒\n', T);
fprintf('  样本数 N = %d\n\n', N);

%% 2. 理论值
theoretical_RMS = A / sqrt(2);
theoretical_Power = theoretical_RMS^2;  % Pa²

fprintf('理论值:\n');
fprintf('  RMS = %.4f Pa\n', theoretical_RMS);
fprintf('  功率 = %.4f Pa²\n\n', theoretical_Power);

%% 3. 使用 computeBandEnergy 计算（1Hz带宽）
fprintf('【方法1】computeBandEnergy（1Hz带宽）\n');
fprintf('----------------------------------------------\n');

% 定义多个1Hz带宽的频带，围绕100Hz
test_freqs = [95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105];
numTestBands = length(test_freqs);

FreqBands_1Hz = zeros(numTestBands, 2);
for i = 1:numTestBands
    FreqBands_1Hz(i, :) = [test_freqs(i), test_freqs(i) + 1];
end

bandPower_1Hz = computeBandEnergy(signal, Fs, FreqBands_1Hz);

fprintf('中心频率 | 频带范围 | 功率(Pa²)\n');
fprintf('---------|----------|----------\n');
for i = 1:numTestBands
    fc = test_freqs(i) + 0.5;
    fprintf('%6.1f Hz | %3d-%3d Hz | %.6f\n', ...
            fc, test_freqs(i), test_freqs(i)+1, bandPower_1Hz(i));
end
fprintf('\n');

%% 4. 使用 pwelch 计算功率谱密度
fprintf('【方法2】pwelch 功率谱密度\n');
fprintf('----------------------------------------------\n');

% pwelch 参数设置
window_length = min(4096, N);  % 窗长度
window = hanning(window_length);
noverlap = floor(window_length / 2);
nfft = max(8192, 2^nextpow2(N));

[Pxx, f_psd] = pwelch(signal, window, noverlap, nfft, Fs);

fprintf('pwelch 参数:\n');
fprintf('  窗长度 = %d\n', window_length);
fprintf('  重叠 = %d\n', noverlap);
fprintf('  NFFT = %d\n', nfft);
fprintf('  频率分辨率 = %.4f Hz\n\n', Fs/nfft);

% 提取对应频率的PSD值
fprintf('中心频率 | PSD(Pa²/Hz) | 功率(PSD×1Hz)\n');
fprintf('---------|-------------|---------------\n');
PSD_at_freqs = zeros(numTestBands, 1);
Power_from_PSD = zeros(numTestBands, 1);

for i = 1:numTestBands
    fc = test_freqs(i) + 0.5;
    [~, idx] = min(abs(f_psd - fc));
    PSD_at_freqs(i) = Pxx(idx);
    Power_from_PSD(i) = PSD_at_freqs(i) * 1;  % PSD × 1Hz带宽
    fprintf('%6.1f Hz | %.6f    | %.6f\n', fc, PSD_at_freqs(i), Power_from_PSD(i));
end
fprintf('\n');

%% 5. 对比结果
fprintf('【对比分析】\n');
fprintf('----------------------------------------------\n');
fprintf('中心频率 | computeBandEnergy | pwelch×1Hz | 差异(dB)\n');
fprintf('---------|-------------------|------------|---------\n');

for i = 1:numTestBands
    fc = test_freqs(i) + 0.5;
    diff_dB = 10 * log10(bandPower_1Hz(i) / Power_from_PSD(i));
    fprintf('%6.1f Hz | %.6f Pa²      | %.6f Pa² | %+6.2f\n', ...
            fc, bandPower_1Hz(i), Power_from_PSD(i), diff_dB);
end

% 计算100Hz附近的平均差异
idx_100Hz = find(test_freqs == 100);
if ~isempty(idx_100Hz)
    diff_at_100Hz = 10 * log10(bandPower_1Hz(idx_100Hz) / Power_from_PSD(idx_100Hz));
    fprintf('\n在100Hz处的差异: %.2f dB\n', diff_at_100Hz);
end

%% 6. 绘制对比图
fprintf('\n【绘制对比图】\n');
fprintf('----------------------------------------------\n');

figure('Position', [100, 100, 1400, 800]);

% 子图1：pwelch 完整频谱
subplot(3, 1, 1);
plot(f_psd, 10*log10(Pxx), 'b-', 'LineWidth', 1);
hold on;
% 标记测试频率
for i = 1:numTestBands
    fc = test_freqs(i) + 0.5;
    [~, idx] = min(abs(f_psd - fc));
    plot(f_psd(idx), 10*log10(Pxx(idx)), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
end
grid on;
xlabel('频率 (Hz)');
ylabel('PSD (dB re 1 Pa²/Hz)');
title('pwelch 功率谱密度（完整频谱）');
xlim([0, 200]);

% 子图2：放大100Hz附近
subplot(3, 1, 2);
zoom_range = (f_psd >= 90) & (f_psd <= 110);
plot(f_psd(zoom_range), 10*log10(Pxx(zoom_range)), 'b-', 'LineWidth', 2);
hold on;

% 叠加 computeBandEnergy 结果（转为等效PSD）
fc_centers = test_freqs + 0.5;
equivalent_PSD_from_computeBandEnergy = bandPower_1Hz;  % 1Hz带宽
stem(fc_centers, 10*log10(equivalent_PSD_from_computeBandEnergy), ...
     'r', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8, ...
     'DisplayName', 'computeBandEnergy (1Hz带宽)');

grid on;
xlabel('频率 (Hz)');
ylabel('功率/PSD (dB)');
title('100Hz附近放大对比');
legend('pwelch PSD', 'computeBandEnergy (等效PSD)');
xlim([90, 110]);

% 子图3：差异曲线
subplot(3, 1, 3);
diff_dB_all = 10 * log10(bandPower_1Hz ./ Power_from_PSD);
stem(fc_centers, diff_dB_all, 'k', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8);
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('频率 (Hz)');
ylabel('差异 (dB)');
title('computeBandEnergy vs pwelch 差异（正值表示computeBandEnergy更大）');
xlim([90, 110]);

fprintf('图形已生成\n\n');

%% 7. 分析差异原因
fprintf('【差异分析】\n');
fprintf('----------------------------------------------\n');

mean_diff = mean(abs(diff_dB_all));
max_diff = max(abs(diff_dB_all));

fprintf('平均绝对差异: %.2f dB\n', mean_diff);
fprintf('最大绝对差异: %.2f dB\n\n', max_diff);

fprintf('可能的差异原因:\n');
fprintf('1. 窗函数:\n');
fprintf('   - computeBandEnergy: 使用矩形窗（无窗）\n');
fprintf('   - pwelch: 使用Hanning窗\n');
fprintf('   → 窗函数会影响功率估计\n\n');

fprintf('2. 频率分辨率:\n');
fprintf('   - computeBandEnergy: 频率分辨率 = %.4f Hz\n', Fs/2^nextpow2(N));
fprintf('   - pwelch: 频率分辨率 = %.4f Hz\n', Fs/nfft);
fprintf('   → 频率分辨率不同导致采样点不同\n\n');

fprintf('3. 分段平均:\n');
fprintf('   - computeBandEnergy: 单次FFT\n');
fprintf('   - pwelch: 多段平均（减少方差）\n');
fprintf('   → pwelch更平滑，噪声更小\n\n');

fprintf('4. 归一化方式:\n');
fprintf('   - computeBandEnergy: 除以样本数N\n');
fprintf('   - pwelch: 有窗函数的能量归一化\n');
fprintf('   → 归一化系数可能不同\n\n');

%% 8. 结论
fprintf('========== 结论 ==========\n');
if mean_diff < 1
    fprintf('✓ 当频带宽度为1Hz时，computeBandEnergy 与 pwelch 结果基本一致\n');
    fprintf('  （平均差异 < 1 dB）\n');
elseif mean_diff < 3
    fprintf('○ computeBandEnergy 与 pwelch 结果接近\n');
    fprintf('  （平均差异 < 3 dB，差异主要由窗函数引起）\n');
else
    fprintf('✗ computeBandEnergy 与 pwelch 存在明显差异\n');
    fprintf('  （平均差异 = %.2f dB）\n', mean_diff);
end

fprintf('\n注意：\n');
fprintf('- computeBandEnergy 适合宽频带（如1/3倍频程）\n');
fprintf('- pwelch 更适合精细频谱分析\n');
fprintf('- 1Hz带宽是一个特殊情况，两者应该接近但不完全相同\n');
fprintf('===============================\n');
