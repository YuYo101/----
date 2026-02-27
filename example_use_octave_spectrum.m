%% 1/3倍频程谱计算示例
% 演示如何使用 compute1_3OctaveSpectrum 函数
% 作者: Claude Code
% 日期: 2026-02-27

clear; clc; close all;

%% 示例1: 基本使用 - 读取音频文件并计算
fprintf('========== 示例1: 基本使用 ==========\n');

% 读取音频文件
audioFile = '10ch_10s_hp30_lp20000_20260209_road_shuinilu_speed_40_03.wav';
[audioData, Fs] = audioread(audioFile);

fprintf('读取音频文件: %s\n', audioFile);
fprintf('  采样频率: %d Hz\n', Fs);
fprintf('  通道数: %d\n', size(audioData, 2));
fprintf('  时长: %.2f 秒\n\n', size(audioData, 1) / Fs);

% 计算1/3倍频程谱
[fc, bandSPL_dB, bandRMS, FreqBands] = compute1_3OctaveSpectrum(audioData, Fs);

% 显示结果（通道1）
fprintf('通道1的1/3倍频程谱:\n');
fprintf('中心频率(Hz) | 下限(Hz) | 上限(Hz) | RMS声压(Pa) | 声压级(dB SPL)\n');
fprintf('-------------|----------|----------|-------------|---------------\n');
for i = 1:length(fc)
    fprintf('%8.1f Hz  | %7.2f  | %7.2f  | %11.4e  | %12.2f\n', ...
            fc(i), FreqBands(i,1), FreqBands(i,2), bandRMS(i,1), bandSPL_dB(i,1));
end

%% 示例2: 自定义频率范围
fprintf('\n========== 示例2: 自定义频率范围 ==========\n');

% 只计算低频部分 (20 Hz - 200 Hz)
[fc_low, SPL_low] = compute1_3OctaveSpectrum(audioData, Fs, ...
                                              'FreqRange', [20, 200]);

fprintf('低频范围 (20-200 Hz):\n');
fprintf('中心频率 | 通道1 SPL | 通道2 SPL\n');
fprintf('---------|-----------|----------\n');
for i = 1:length(fc_low)
    fprintf('%6.1f Hz | %8.2f dB | %8.2f dB\n', ...
            fc_low(i), SPL_low(i,1), SPL_low(i,2));
end

%% 示例3: 绘制对比图
fprintf('\n========== 示例3: 绘制多通道对比 ==========\n');

figure('Position', [100, 100, 1200, 500]);

% 选择前4个通道绘制
numPlotChannels = min(4, size(audioData, 2));
colors = lines(numPlotChannels);

hold on;
for ch = 1:numPlotChannels
    semilogx(fc, bandSPL_dB(:, ch), '-o', ...
             'LineWidth', 2, 'MarkerSize', 6, ...
             'Color', colors(ch,:), ...
             'DisplayName', sprintf('通道 %d', ch));
end

grid on;
xlabel('中心频率 (Hz)', 'FontSize', 12);
ylabel('声压级 (dB SPL)', 'FontSize', 12);
title('1/3倍频程谱对比', 'FontSize', 14);
legend('Location', 'best');
xlim([fc(1)*0.9, fc(end)*1.1]);

fprintf('图形已生成\n');

%% 示例4: 计算总声压级
fprintf('\n========== 示例4: 计算总声压级 ==========\n');

% 计算各频带能量
N = size(audioData, 1);
bandEnergy = bandRMS.^2 * N;  % 从RMS反推能量

% 计算总能量和总声压级
for ch = 1:size(audioData, 2)
    totalEnergy = sum(bandEnergy(:, ch));
    totalRMS = sqrt(totalEnergy / N);
    totalSPL = 20 * log10(totalRMS / 2e-5);

    fprintf('通道 %2d: 总RMS = %.4f Pa, 总SPL = %.2f dB SPL\n', ...
            ch, totalRMS, totalSPL);
end

%% 示例5: 导出为CSV
fprintf('\n========== 示例5: 导出为CSV ==========\n');

csvFile = 'example_octave_spectrum.csv';
fid = fopen(csvFile, 'w');

% 写入表头
fprintf(fid, '中心频率(Hz),下限(Hz),上限(Hz)');
for ch = 1:size(bandSPL_dB, 2)
    fprintf(fid, ',通道%d_SPL(dB)', ch);
end
fprintf(fid, '\n');

% 写入数据
for i = 1:length(fc)
    fprintf(fid, '%.1f,%.2f,%.2f', fc(i), FreqBands(i,1), FreqBands(i,2));
    for ch = 1:size(bandSPL_dB, 2)
        fprintf(fid, ',%.4f', bandSPL_dB(i, ch));
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('结果已保存到: %s\n', csvFile);

%% 示例6: 找出峰值频率
fprintf('\n========== 示例6: 找出峰值频率 ==========\n');

for ch = 1:min(3, size(bandSPL_dB, 2))
    [maxSPL, maxIdx] = max(bandSPL_dB(:, ch));
    fprintf('通道 %d: 峰值频率 = %.1f Hz, 声压级 = %.2f dB SPL\n', ...
            ch, fc(maxIdx), maxSPL);
end

fprintf('\n========== 所有示例完成 ==========\n');
