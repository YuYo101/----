%% 对比 computeBandEnergy 与 computeBandPower_pwelch
% 展示两种方法的差异及各自的优缺点
clear; clc; close all;

fprintf('========== 两种频带功率计算方法对比 ==========\n\n');

%% 1. 生成测试信号
Fs = 8000;  % 8 kHz
T = 2;      % 2 秒
t = (0:1/Fs:T-1/Fs)';
N = length(t);

% 多频率成分 + 噪声
signal = 2.0*sin(2*pi*40*t) + ...     % 40 Hz, 2 Pa
         1.5*sin(2*pi*100*t) + ...    % 100 Hz, 1.5 Pa
         1.0*sin(2*pi*1000*t) + ...   % 1000 Hz, 1 Pa
         0.05*randn(N, 1);            % 白噪声

fprintf('测试信号:\n');
fprintf('  40 Hz: 2.0 Pa (RMS=%.3f Pa)\n', 2/sqrt(2));
fprintf('  100 Hz: 1.5 Pa (RMS=%.3f Pa)\n', 1.5/sqrt(2));
fprintf('  1000 Hz: 1.0 Pa (RMS=%.3f Pa)\n', 1/sqrt(2));
fprintf('  + 白噪声 (σ=0.05)\n');
fprintf('  采样频率: %d Hz\n', Fs);
fprintf('  样本数: %d\n\n', N);

%% 2. 定义1/3倍频程频带
fc_octave = [20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, ...
             250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000];
bandwidth_ratio = 2^(1/6);
fl = fc_octave / bandwidth_ratio;
fu = fc_octave * bandwidth_ratio;
FreqBands = [fl', fu'];

fprintf('频带定义: 1/3倍频程（20 Hz ~ 2000 Hz）\n');
fprintf('频带数: %d\n\n', length(fc_octave));

%% 3. 方法1：computeBandEnergy（基于FFT）
fprintf('【方法1】computeBandEnergy (单次FFT + 矩形窗)\n');
fprintf('----------------------------------------------\n');
tic;
bandPower_FFT = computeBandEnergy(signal, Fs, FreqBands);
time_FFT = toc;

bandRMS_FFT = sqrt(bandPower_FFT);
bandSPL_FFT = 20 * log10(bandRMS_FFT / 2e-5);

fprintf('计算完成，耗时: %.4f 秒\n\n', time_FFT);

%% 4. 方法2：computeBandPower_pwelch（基于pwelch）
fprintf('【方法2】computeBandPower_pwelch (pwelch + Hanning窗)\n');
fprintf('----------------------------------------------\n');
tic;
bandPower_pwelch = computeBandPower_pwelch(signal, Fs, FreqBands);
time_pwelch = toc;

bandRMS_pwelch = sqrt(bandPower_pwelch);
bandSPL_pwelch = 20 * log10(bandRMS_pwelch / 2e-5);

fprintf('计算完成，耗时: %.4f 秒\n\n', time_pwelch);

%% 5. 详细对比
fprintf('【详细对比】\n');
fprintf('----------------------------------------------\n');
fprintf('中心频率 | 方法1功率  | 方法2功率  | 方法1 SPL | 方法2 SPL | 差异(dB)\n');
fprintf('---------|-----------|-----------|-----------|-----------|----------\n');

for i = 1:length(fc_octave)
    diff_dB = bandSPL_FFT(i) - bandSPL_pwelch(i);
    fprintf('%6.1f Hz | %.4e | %.4e | %7.2f dB | %7.2f dB | %+6.2f\n', ...
            fc_octave(i), bandPower_FFT(i), bandPower_pwelch(i), ...
            bandSPL_FFT(i), bandSPL_pwelch(i), diff_dB);
end

% 统计差异
diff_dB_all = bandSPL_FFT - bandSPL_pwelch;
mean_diff = mean(abs(diff_dB_all));
max_diff = max(abs(diff_dB_all));
std_diff = std(diff_dB_all);

fprintf('\n统计:\n');
fprintf('  平均绝对差异: %.2f dB\n', mean_diff);
fprintf('  最大绝对差异: %.2f dB\n', max_diff);
fprintf('  标准差: %.2f dB\n\n', std_diff);

%% 6. 重点频率对比
fprintf('【重点频率对比】（包含信号成分）\n');
fprintf('----------------------------------------------\n');
key_freqs = [40, 100, 1000];
fprintf('频率    | 理论RMS | 方法1 RMS | 方法2 RMS | 方法1误差 | 方法2误差\n');
fprintf('--------|---------|-----------|-----------|-----------|----------\n');

for f = key_freqs
    [~, idx] = min(abs(fc_octave - f));
    theoretical_RMS = f/40 * 2/sqrt(2);  % 简化：按比例估算
    if f == 40
        theoretical_RMS = 2/sqrt(2);
    elseif f == 100
        theoretical_RMS = 1.5/sqrt(2);
    elseif f == 1000
        theoretical_RMS = 1/sqrt(2);
    end

    err1 = abs(bandRMS_FFT(idx) - theoretical_RMS) / theoretical_RMS * 100;
    err2 = abs(bandRMS_pwelch(idx) - theoretical_RMS) / theoretical_RMS * 100;

    fprintf('%4d Hz | %.4f  | %.4f    | %.4f    | %6.2f%%   | %6.2f%%\n', ...
            f, theoretical_RMS, bandRMS_FFT(idx), bandRMS_pwelch(idx), err1, err2);
end
fprintf('\n');

%% 7. 绘制对比图
fprintf('【绘制对比图】\n');
fprintf('----------------------------------------------\n');

figure('Position', [100, 100, 1400, 900]);

% 子图1：声压级对比
subplot(3, 2, 1);
semilogx(fc_octave, bandSPL_FFT, 'b-o', 'LineWidth', 2, 'MarkerSize', 6, ...
         'DisplayName', 'computeBandEnergy (FFT)');
hold on;
semilogx(fc_octave, bandSPL_pwelch, 'r--s', 'LineWidth', 2, 'MarkerSize', 6, ...
         'DisplayName', 'computeBandPower\_pwelch');
grid on;
xlabel('中心频率 (Hz)');
ylabel('声压级 (dB SPL)');
title('1/3倍频程谱对比');
legend('Location', 'best');
xlim([fc_octave(1)*0.9, fc_octave(end)*1.1]);

% 子图2：差异曲线
subplot(3, 2, 2);
semilogx(fc_octave, diff_dB_all, 'k-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
yline(mean_diff, 'g--', sprintf('平均: %.2f dB', mean_diff));
yline(-mean_diff, 'g--', sprintf('平均: -%.2f dB', mean_diff));
grid on;
xlabel('中心频率 (Hz)');
ylabel('差异 (dB)');
title('声压级差异 (FFT - pwelch)');
xlim([fc_octave(1)*0.9, fc_octave(end)*1.1]);

% 子图3：功率对比（线性）
subplot(3, 2, 3);
bar_data = [bandPower_FFT, bandPower_pwelch];
bar(1:length(fc_octave), bar_data);
set(gca, 'XTick', 1:2:length(fc_octave));
set(gca, 'XTickLabel', fc_octave(1:2:end));
xlabel('中心频率 (Hz)');
ylabel('功率 (Pa²)');
title('频带功率对比（线性尺度）');
legend('FFT', 'pwelch', 'Location', 'best');
grid on;

% 子图4：相对误差
subplot(3, 2, 4);
rel_error = abs(bandPower_FFT - bandPower_pwelch) ./ bandPower_pwelch * 100;
semilogx(fc_octave, rel_error, 'm-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
yline(10, 'r--', '10%');
yline(5, 'g--', '5%');
grid on;
xlabel('中心频率 (Hz)');
ylabel('相对误差 (%)');
title('功率相对误差');
xlim([fc_octave(1)*0.9, fc_octave(end)*1.1]);

% 子图5：噪声底对比（远离信号频率）
subplot(3, 2, 5);
noise_bands = [1:3, 7:9, 11:14, 16:18];  % 远离40, 100, 1000 Hz的频带
semilogx(fc_octave(noise_bands), bandSPL_FFT(noise_bands), 'b-o', ...
         'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'FFT');
hold on;
semilogx(fc_octave(noise_bands), bandSPL_pwelch(noise_bands), 'r--s', ...
         'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'pwelch');
grid on;
xlabel('中心频率 (Hz)');
ylabel('声压级 (dB SPL)');
title('噪声底对比（远离信号频率）');
legend('Location', 'best');

% 子图6：计算时间对比
subplot(3, 2, 6);
bar_times = [time_FFT, time_pwelch];
bar(bar_times);
set(gca, 'XTickLabel', {'FFT', 'pwelch'});
ylabel('时间 (秒)');
title(sprintf('计算时间对比（FFT快%.1f倍）', time_pwelch/time_FFT));
grid on;
for i = 1:2
    text(i, bar_times(i), sprintf('%.4fs', bar_times(i)), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end

fprintf('图形已生成\n\n');

%% 8. 方法对比总结
fprintf('========== 方法对比总结 ==========\n');
fprintf('\n【computeBandEnergy (FFT方法)】\n');
fprintf('优点:\n');
fprintf('  ✓ 速度快（%.4f秒，快%.1f倍）\n', time_FFT, time_pwelch/time_FFT);
fprintf('  ✓ 实时性好\n');
fprintf('  ✓ 适合宽频带（如1/3倍频程）\n');
fprintf('缺点:\n');
fprintf('  ✗ 频率泄漏较大（矩形窗）\n');
fprintf('  ✗ 噪声方差较大（单次估计）\n');
fprintf('  ✗ 噪声底估计偏高\n');

fprintf('\n【computeBandPower_pwelch (Welch方法)】\n');
fprintf('优点:\n');
fprintf('  ✓ 频率泄漏小（Hanning窗）\n');
fprintf('  ✓ 噪声方差小（分段平均）\n');
fprintf('  ✓ 精度高（尤其在噪声区域）\n');
fprintf('缺点:\n');
fprintf('  ✗ 计算稍慢（%.4f秒）\n', time_pwelch);
fprintf('  ✗ 窗函数损失部分频率分辨率\n');

fprintf('\n【使用建议】\n');
if mean_diff < 2
    fprintf('✓ 在主要频率成分处，两种方法结果接近（差异<2dB）\n');
    fprintf('  → 宽频带应用可以使用computeBandEnergy（更快）\n');
else
    fprintf('○ 两种方法存在一定差异（平均%.2fdB）\n', mean_diff);
end
fprintf('  → 需要高精度/低噪声: 使用computeBandPower_pwelch\n');
fprintf('  → 需要实时性/高速度: 使用computeBandEnergy\n');
fprintf('  → 噪声分析/精细谱: 使用computeBandPower_pwelch\n');
fprintf('  → 1/3倍频程分析: 两者均可，FFT更快\n');

fprintf('\n================================\n');
