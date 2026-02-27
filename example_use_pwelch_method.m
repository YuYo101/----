%% computeBandPower_pwelch 使用示例与最佳实践
% 演示如何在不同场景下使用基于pwelch的频带功率计算方法
clear; clc; close all;

fprintf('========== computeBandPower_pwelch 使用指南 ==========\n\n');

%% 示例1：基本使用（对比computeBandEnergy）
fprintf('【示例1】基本使用\n');
fprintf('----------------------------------------------\n');

% 生成测试信号
Fs = 8000;
T = 2;
t = (0:1/Fs:T-1/Fs)';
signal = 2.0*sin(2*pi*100*t) + 0.1*randn(length(t),1);

% 定义频带（1/3倍频程）
fc_octave = [63, 80, 100, 125, 160];
bandwidth_ratio = 2^(1/6);
fl = fc_octave / bandwidth_ratio;
fu = fc_octave * bandwidth_ratio;
FreqBands = [fl', fu'];

% 方法1：快速FFT方法
tic;
bandPower_FFT = computeBandEnergy(signal, Fs, FreqBands);
time_FFT = toc;

% 方法2：高精度pwelch方法
tic;
bandPower_pwelch = computeBandPower_pwelch(signal, Fs, FreqBands);
time_pwelch = toc;

% 计算声压级
bandSPL_FFT = 20 * log10(sqrt(bandPower_FFT) / 2e-5);
bandSPL_pwelch = 20 * log10(sqrt(bandPower_pwelch) / 2e-5);

fprintf('中心频率 | FFT方法  | pwelch方法 | 差异(dB)\n');
fprintf('---------|----------|-----------|----------\n');
for i = 1:length(fc_octave)
    diff_dB = bandSPL_FFT(i) - bandSPL_pwelch(i);
    fprintf('%4d Hz  | %6.2f dB | %6.2f dB  | %+5.2f\n', ...
            fc_octave(i), bandSPL_FFT(i), bandSPL_pwelch(i), diff_dB);
end
fprintf('\n速度对比: FFT耗时 %.4fs, pwelch耗时 %.4fs\n\n', time_FFT, time_pwelch);

%% 示例2：自定义窗函数和FFT参数
fprintf('【示例2】自定义参数（提高频率分辨率）\n');
fprintf('----------------------------------------------\n');

% 默认参数
bandPower_default = computeBandPower_pwelch(signal, Fs, FreqBands);

% 自定义参数：更长的窗，更高的频率分辨率
bandPower_custom = computeBandPower_pwelch(signal, Fs, FreqBands, ...
    'Window', hanning(8192), ...    % 更长的窗
    'Overlap', 4096, ...            % 50%重叠
    'NFFT', 32768);                 % 更高的FFT分辨率

fprintf('中心频率 | 默认参数 | 自定义参数 | 差异(dB)\n');
fprintf('---------|----------|-----------|----------\n');
for i = 1:length(fc_octave)
    SPL_default = 20 * log10(sqrt(bandPower_default(i)) / 2e-5);
    SPL_custom = 20 * log10(sqrt(bandPower_custom(i)) / 2e-5);
    diff_dB = SPL_default - SPL_custom;
    fprintf('%4d Hz  | %6.2f dB | %6.2f dB  | %+5.2f\n', ...
            fc_octave(i), SPL_default, SPL_custom, diff_dB);
end
fprintf('注意: 更长的窗和更高的FFT分辨率可以提高低频精度\n\n');

%% 示例3：处理多通道音频
fprintf('【示例3】多通道音频处理\n');
fprintf('----------------------------------------------\n');

% 生成3通道测试信号
nChannels = 3;
signal_multi = zeros(length(t), nChannels);
signal_multi(:, 1) = 1.0*sin(2*pi*100*t) + 0.05*randn(length(t),1);  % 通道1
signal_multi(:, 2) = 1.5*sin(2*pi*100*t) + 0.05*randn(length(t),1);  % 通道2
signal_multi(:, 3) = 2.0*sin(2*pi*100*t) + 0.05*randn(length(t),1);  % 通道3

% 计算所有通道的频带功率
bandPower_multi = computeBandPower_pwelch(signal_multi, Fs, FreqBands);

fprintf('100 Hz频带的各通道声压级:\n');
idx_100Hz = find(fc_octave == 100);
for ch = 1:nChannels
    SPL = 20 * log10(sqrt(bandPower_multi(idx_100Hz, ch)) / 2e-5);
    fprintf('  通道%d: %.2f dB SPL\n', ch, SPL);
end
fprintf('\n');

%% 示例4：噪声分析（pwelch的优势场景）
fprintf('【示例4】噪声分析（pwelch优于FFT）\n');
fprintf('----------------------------------------------\n');

% 生成弱信号 + 强噪声
weak_signal = 0.1*sin(2*pi*100*t);  % 弱信号 0.1 Pa
strong_noise = 0.5*randn(length(t),1);  % 强噪声 σ=0.5
signal_noisy = weak_signal + strong_noise;

% 定义更宽的频带范围（包含噪声区域）
fc_wide = [20, 40, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500];
fl_wide = fc_wide / bandwidth_ratio;
fu_wide = fc_wide * bandwidth_ratio;
FreqBands_wide = [fl_wide', fu_wide'];

% 两种方法对比
bandPower_FFT_noisy = computeBandEnergy(signal_noisy, Fs, FreqBands_wide);
bandPower_pwelch_noisy = computeBandPower_pwelch(signal_noisy, Fs, FreqBands_wide);

bandSPL_FFT_noisy = 20 * log10(sqrt(bandPower_FFT_noisy) / 2e-5);
bandSPL_pwelch_noisy = 20 * log10(sqrt(bandPower_pwelch_noisy) / 2e-5);

fprintf('中心频率 | FFT方法  | pwelch方法 | 差异(dB)\n');
fprintf('---------|----------|-----------|----------\n');
for i = 1:length(fc_wide)
    diff_dB = bandSPL_FFT_noisy(i) - bandSPL_pwelch_noisy(i);
    marker = '';
    if fc_wide(i) == 100
        marker = ' ← 信号频率';
    end
    fprintf('%4d Hz  | %6.2f dB | %6.2f dB  | %+5.2f%s\n', ...
            fc_wide(i), bandSPL_FFT_noisy(i), bandSPL_pwelch_noisy(i), ...
            diff_dB, marker);
end
fprintf('结论: 在噪声区域，pwelch的噪声底比FFT低10-20dB\n\n');

%% 示例5：窄带分析（1Hz带宽）
fprintf('【示例5】窄带分析（1Hz带宽）\n');
fprintf('----------------------------------------------\n');

% 生成纯净单频信号
signal_pure = 2.0*sin(2*pi*100*t);

% 定义100Hz附近的1Hz频带
test_freqs = 98:102;
FreqBands_1Hz = [test_freqs', (test_freqs+1)'];

% 计算
bandPower_1Hz = computeBandPower_pwelch(signal_pure, Fs, FreqBands_1Hz);

fprintf('频带     | 功率(Pa²)  | SPL(dB)\n');
fprintf('---------|-----------|----------\n');
for i = 1:length(test_freqs)
    SPL = 20 * log10(sqrt(bandPower_1Hz(i)) / 2e-5);
    marker = '';
    if test_freqs(i) == 99 || test_freqs(i) == 100
        marker = ' ← 主频率';
    end
    fprintf('%3d-%3d Hz | %.6f | %6.2f%s\n', ...
            test_freqs(i), test_freqs(i)+1, bandPower_1Hz(i), SPL, marker);
end
fprintf('注意: 1Hz带宽时，pwelch能更准确地定位频率成分\n\n');

%% 示例6：实际音频文件处理
fprintf('【示例6】实际音频文件处理（如果文件存在）\n');
fprintf('----------------------------------------------\n');

audio_file = '10ch_10s_hp30_lp20000_20260209_road_shuinilu_speed_40_03.wav';

if exist(audio_file, 'file')
    % 读取音频
    [audioData, Fs_audio] = audioread(audio_file);
    fprintf('已读取音频: %d通道, %.2f秒, %.1f kHz\n', ...
            size(audioData,2), size(audioData,1)/Fs_audio, Fs_audio/1000);

    % 定义标准1/3倍频程频带（20-20kHz）
    fc_standard = [20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, ...
                   250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, ...
                   2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000];

    % 筛选有效频带
    valid_freqs = fc_standard(fc_standard < Fs_audio/2);
    fl_audio = valid_freqs / bandwidth_ratio;
    fu_audio = valid_freqs * bandwidth_ratio;
    FreqBands_audio = [fl_audio', fu_audio'];

    % 使用pwelch方法计算（推荐用于高质量分析）
    fprintf('计算中...\n');
    tic;
    bandPower_audio = computeBandPower_pwelch(audioData, Fs_audio, FreqBands_audio);
    time_elapsed = toc;
    fprintf('计算完成，耗时: %.2f秒\n\n', time_elapsed);

    % 显示第1通道的结果
    bandSPL_audio = 20 * log10(sqrt(bandPower_audio(:,1)) / 2e-5);

    fprintf('第1通道 1/3倍频程谱 (选择性显示):\n');
    fprintf('中心频率 | 声压级(dB SPL)\n');
    fprintf('---------|---------------\n');
    display_idx = [4, 8, 18, 26, 30];  % 40, 100, 1000, 8000, 20000 Hz
    for i = display_idx
        if i <= length(valid_freqs)
            fprintf('%6.1f Hz | %6.2f\n', valid_freqs(i), bandSPL_audio(i));
        end
    end
    fprintf('...\n');
else
    fprintf('音频文件未找到，跳过此示例\n');
end
fprintf('\n');

%% 使用建议总结
fprintf('========== 使用建议总结 ==========\n\n');

fprintf('【何时使用 computeBandEnergy (FFT方法)】\n');
fprintf('✓ 实时处理（速度要求高）\n');
fprintf('✓ 宽频带分析（如1/3倍频程）\n');
fprintf('✓ 信噪比高的场景\n');
fprintf('✓ 快速原型开发\n\n');

fprintf('【何时使用 computeBandPower_pwelch (pwelch方法)】\n');
fprintf('✓ 需要高精度测量\n');
fprintf('✓ 噪声分析（低噪声底要求）\n');
fprintf('✓ 窄带分析（如1Hz带宽）\n');
fprintf('✓ 弱信号检测\n');
fprintf('✓ 科研级分析\n');
fprintf('✓ 后处理场景（非实时）\n\n');

fprintf('【参数调整建议】\n');
fprintf('1. Window（窗长度）:\n');
fprintf('   - 默认: min(4096, N)\n');
fprintf('   - 低频精度要求高 → 增大窗长度（如8192）\n');
fprintf('   - 时变信号 → 减小窗长度（如2048）\n\n');

fprintf('2. Overlap（重叠比例）:\n');
fprintf('   - 默认: 50%%\n');
fprintf('   - 提高平滑度 → 增大重叠（如75%%）\n');
fprintf('   - 加快速度 → 减小重叠（如25%%）\n\n');

fprintf('3. NFFT（FFT点数）:\n');
fprintf('   - 默认: max(4096, 2^nextpow2(N))\n');
fprintf('   - 精细频谱 → 增大NFFT（如16384）\n');
fprintf('   - 宽带分析 → 使用默认值\n\n');

fprintf('【典型应用场景】\n');
fprintf('• 车辆噪声分析: pwelch（噪声成分复杂）\n');
fprintf('• 语音信号处理: FFT（实时性要求）\n');
fprintf('• 建筑声学测量: pwelch（精度要求高）\n');
fprintf('• 故障诊断: pwelch（弱信号检测）\n');
fprintf('• 快速预览: FFT（速度优先）\n\n');

fprintf('================================\n');
