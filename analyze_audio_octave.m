%% 音频文件1/3倍频程谱分析
% 读取音频文件并计算1/3倍频程谱
% 作者: Claude Code
% 日期: 2026-02-27

clear; clc; close all;

%% 1. 读取音频文件
fprintf('========== 音频文件1/3倍频程谱分析 ==========\n\n');
fprintf('【步骤1】读取音频文件\n');
fprintf('----------------------------------------------\n');

audioFile = '10ch_10s_hp30_lp20000_20260209_road_shuinilu_speed_40_03.wav';

% 检查文件是否存在
if ~exist(audioFile, 'file')
    error('音频文件不存在: %s', audioFile);
end

% 读取音频信息
info = audioinfo(audioFile);
fprintf('音频文件信息:\n');
fprintf('  文件名: %s\n', audioFile);
fprintf('  采样频率: %d Hz\n', info.SampleRate);
fprintf('  通道数: %d\n', info.NumChannels);
fprintf('  时长: %.2f 秒\n', info.Duration);
fprintf('  样本数: %d\n', info.TotalSamples);
fprintf('  位深度: %d bit\n', info.BitsPerSample);

% 读取音频数据
fprintf('\n正在读取音频数据...\n');
[audioData, Fs] = audioread(audioFile);
[N, nChannels] = size(audioData);
fprintf('读取完成！数据尺寸: %d x %d\n\n', N, nChannels);

%% 2. 定义1/3倍频程频带
fprintf('【步骤2】定义1/3倍频程频带\n');
fprintf('----------------------------------------------\n');

% 标准1/3倍频程中心频率（IEC 61260标准）
fc_standard = [20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, ...
               250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, ...
               2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000];

% 根据Nyquist频率筛选有效频带
fc_valid = fc_standard(fc_standard <= Fs/2);

% 计算上下限频率
bandwidth_ratio = 2^(1/6);  % 2^(1/6) ≈ 1.122
fl = fc_valid / bandwidth_ratio;
fu = fc_valid * bandwidth_ratio;

% 构造频带矩阵
FreqBands = [fl(:), fu(:)];
numBands = length(fc_valid);

fprintf('采样频率: %d Hz\n', Fs);
fprintf('Nyquist频率: %d Hz\n', Fs/2);
fprintf('有效频带数量: %d (%.1f Hz ~ %.0f Hz)\n\n', numBands, fc_valid(1), fc_valid(end));

%% 3. 计算1/3倍频程谱（每个通道）
fprintf('【步骤3】计算1/3倍频程谱\n');
fprintf('----------------------------------------------\n');

tic;
bandEnergy = computeBandEnergy(audioData, Fs, FreqBands);
time_elapsed = toc;

fprintf('计算完成，耗时: %.4f 秒\n', time_elapsed);
fprintf('结果尺寸: %d x %d (频带数 x 通道数)\n\n', size(bandEnergy, 1), size(bandEnergy, 2));

%% 4. 转换为声压级（SPL, dB）
fprintf('【步骤4】转换为声压级（SPL, dB）\n');
fprintf('----------------------------------------------\n');

% 计算各频带的RMS声压值
% RMS = sqrt(能量 / 样本数)
bandRMS = sqrt(bandEnergy / N);  % 均方根值（假设音频单位是Pa）

% 声压参考值（空气中的标准参考声压）
p_ref = 2e-5;  % 20 μPa = 2×10^-5 Pa

% 计算声压级（SPL）
% SPL = 20 * log10(p_rms / p_ref)
bandSPL_dB = 20 * log10(bandRMS / p_ref + eps);

fprintf('转换完成\n');
fprintf('参考声压: %.2e Pa (20 μPa)\n', p_ref);
fprintf('SPL范围: %.2f ~ %.2f dB SPL\n\n', min(bandSPL_dB(:)), max(bandSPL_dB(:)));

%% 5. 显示各通道的1/3倍频程谱
fprintf('【步骤5】显示1/3倍频程谱结果\n');
fprintf('----------------------------------------------\n');

for ch = 1:nChannels
    fprintf('\n通道 %d 的1/3倍频程谱:\n', ch);
    fprintf('中心频率(Hz) |  RMS声压(Pa)  | 声压级(dB SPL)\n');
    fprintf('-------------|---------------|---------------\n');

    for i = 1:numBands
        fprintf('%8.1f Hz  | %13.4e | %12.2f\n', ...
                fc_valid(i), bandRMS(i, ch), bandSPL_dB(i, ch));
    end
end

%% 6. 统计信息
fprintf('\n【步骤6】统计信息\n');
fprintf('----------------------------------------------\n');

fprintf('\n各通道总声压级统计:\n');
fprintf('通道 | 总能量(线性) | 总RMS声压(Pa) | 总声压级(dB SPL)\n');
fprintf('-----|--------------|---------------|----------------\n');
for ch = 1:nChannels
    totalEnergy = sum(bandEnergy(:, ch));
    totalRMS = sqrt(totalEnergy / N);
    totalSPL_dB = 20 * log10(totalRMS / p_ref + eps);
    fprintf('  %d  | %12.4e | %13.4e | %13.2f\n', ch, totalEnergy, totalRMS, totalSPL_dB);
end

% 找出主要能量频带（前5个）
fprintf('\n各通道主要能量频带（前5个）:\n');
for ch = 1:nChannels
    fprintf('\n通道 %d:\n', ch);
    [sorted_SPL, idx] = sort(bandSPL_dB(:, ch), 'descend');
    fprintf('  排名 | 中心频率 | 声压级\n');
    fprintf('  -----|----------|-------------\n');
    for i = 1:min(5, numBands)
        fprintf('   %d   | %6.1f Hz | %8.2f dB SPL\n', ...
                i, fc_valid(idx(i)), sorted_SPL(i));
    end
end

%% 7. 绘制频谱图
fprintf('\n【步骤7】绘制频谱图\n');
fprintf('----------------------------------------------\n');

% 如果通道数较多，分两图显示
if nChannels <= 4
    % 少通道：单图显示所有通道
    figure('Position', [100, 100, 1200, 600]);

    colors = lines(nChannels);
    hold on;
    for ch = 1:nChannels
        semilogx(fc_valid, bandSPL_dB(:, ch), '-o', ...
                 'LineWidth', 2, 'MarkerSize', 6, ...
                 'Color', colors(ch,:), ...
                 'DisplayName', sprintf('通道 %d', ch));
    end

    grid on;
    xlabel('中心频率 (Hz)', 'FontSize', 12);
    ylabel('声压级 (dB SPL)', 'FontSize', 12);
    title('1/3倍频程谱 - 所有通道', 'FontSize', 14);
    legend('Location', 'best');
    xlim([min(fc_valid)*0.9, max(fc_valid)*1.1]);

else
    % 多通道：分成多个子图
    numRows = ceil(nChannels / 2);
    figure('Position', [100, 100, 1400, 300*numRows]);

    for ch = 1:nChannels
        subplot(numRows, 2, ch);
        semilogx(fc_valid, bandSPL_dB(:, ch), '-o', ...
                 'LineWidth', 1.5, 'MarkerSize', 5);
        grid on;
        xlabel('中心频率 (Hz)');
        ylabel('声压级 (dB SPL)');
        title(sprintf('通道 %d - 1/3倍频程谱', ch));
        xlim([min(fc_valid)*0.9, max(fc_valid)*1.1]);
    end
end

fprintf('图形已生成\n');

% 绘制全部通道对比热力图（如果通道数>1）
if nChannels > 1
    figure('Position', [100, 150, 1200, 600]);

    imagesc(1:nChannels, 1:numBands, bandSPL_dB);
    colorbar;
    colormap('jet');

    % 设置Y轴刻度为频率
    set(gca, 'YDir', 'normal');
    set(gca, 'YTick', 1:numBands);
    set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%.1f', x), fc_valid, 'UniformOutput', false));

    xlabel('通道号', 'FontSize', 12);
    ylabel('中心频率 (Hz)', 'FontSize', 12);
    title('1/3倍频程谱 - 通道对比热力图 (dB SPL)', 'FontSize', 14);

    % 添加数值标注（可选，如果频带数不太多）
    if numBands <= 20
        for ch = 1:nChannels
            for band = 1:numBands
                text(ch, band, sprintf('%.1f', bandSPL_dB(band, ch)), ...
                     'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', ...
                     'FontSize', 7, ...
                     'Color', 'w');
            end
        end
    end

    fprintf('热力图已生成\n');
end

%% 8. 保存结果
fprintf('\n【步骤8】保存结果\n');
fprintf('----------------------------------------------\n');

% 保存为MAT文件
resultFile = 'octave_spectrum_result.mat';
save(resultFile, 'fc_valid', 'bandEnergy', 'bandRMS', 'bandSPL_dB', ...
     'Fs', 'nChannels', 'audioFile', 'FreqBands', 'p_ref');
fprintf('结果已保存到: %s\n', resultFile);

% 保存为CSV文件（便于导出）
csvFile = 'octave_spectrum_result.csv';
fid = fopen(csvFile, 'w');

% 写入表头
fprintf(fid, '中心频率(Hz),下限(Hz),上限(Hz)');
for ch = 1:nChannels
    fprintf(fid, ',通道%d_声压级(dB_SPL)', ch);
end
fprintf(fid, '\n');

% 写入数据
for i = 1:numBands
    fprintf(fid, '%.1f,%.2f,%.2f', fc_valid(i), fl(i), fu(i));
    for ch = 1:nChannels
        fprintf(fid, ',%.4f', bandSPL_dB(i, ch));
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('结果已保存到: %s\n', csvFile);

fprintf('\n========== 分析完成 ==========\n');
fprintf('所有通道的1/3倍频程谱已计算完成\n');
fprintf('图形和数据文件已生成\n');
fprintf('==============================\n');
