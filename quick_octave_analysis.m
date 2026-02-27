%% 快速计算脚本 - 计算指定音频文件的1/3倍频程谱
% 这是一个简化的计算脚本，可以快速计算任意音频文件的1/3倍频程谱
%
% 使用方法:
%   1. 修改下面的 audioFile 变量为你的音频文件路径
%   2. 运行脚本
%   3. 结果自动保存为 CSV 和 MAT 文件
%
% 作者: Claude Code
% 日期: 2026-02-27

clear; clc; close all;

%% ========== 配置参数（根据需要修改） ==========
audioFile = '10ch_10s_hp30_lp20000_20260209_road_shuinilu_speed_40_03.wav';  % 音频文件路径
freqRange = [20, 20000];      % 频率范围 (Hz)
p_ref = 2e-5;                 % 参考声压 (Pa)，默认20 μPa
outputPrefix = 'octave_';     % 输出文件前缀

%% ========== 1. 读取音频文件 ==========
fprintf('正在读取音频文件...\n');
[audioData, Fs] = audioread(audioFile);
[N, nChannels] = size(audioData);

fprintf('文件: %s\n', audioFile);
fprintf('采样频率: %d Hz\n', Fs);
fprintf('通道数: %d\n', nChannels);
fprintf('时长: %.2f 秒\n', N/Fs);
fprintf('样本数: %d\n\n', N);

%% ========== 2. 计算1/3倍频程谱 ==========
fprintf('正在计算1/3倍频程谱...\n');
tic;
[fc, bandSPL_dB, bandRMS, FreqBands] = compute1_3OctaveSpectrum(audioData, Fs, ...
                                                                 'FreqRange', freqRange, ...
                                                                 'RefPressure', p_ref);
elapsedTime = toc;
fprintf('计算完成！耗时: %.4f 秒\n', elapsedTime);
fprintf('频带数: %d\n', length(fc));
fprintf('SPL范围: %.2f ~ %.2f dB SPL\n\n', min(bandSPL_dB(:)), max(bandSPL_dB(:)));

%% ========== 3. 计算总声压级 ==========
fprintf('各通道总声压级:\n');
fprintf('通道 | 总RMS声压(Pa) | 总声压级(dB SPL)\n');
fprintf('-----|---------------|----------------\n');

bandEnergy = bandRMS.^2 * N;
for ch = 1:nChannels
    totalEnergy = sum(bandEnergy(:, ch));
    totalRMS = sqrt(totalEnergy / N);
    totalSPL = 20 * log10(totalRMS / p_ref);
    fprintf(' %2d  | %13.4e | %13.2f\n', ch, totalRMS, totalSPL);
end
fprintf('\n');

%% ========== 4. 显示主要频带（前5个最高声压级） ==========
fprintf('各通道主要频带（前5个）:\n');
for ch = 1:nChannels
    fprintf('\n通道 %d:\n', ch);
    [sortedSPL, idx] = sort(bandSPL_dB(:, ch), 'descend');
    fprintf('  排名 | 中心频率 | 声压级\n');
    fprintf('  -----|----------|-------------\n');
    for i = 1:min(5, length(fc))
        fprintf('   %d   | %6.1f Hz | %8.2f dB SPL\n', ...
                i, fc(idx(i)), sortedSPL(i));
    end
end
fprintf('\n');

%% ========== 5. 保存结果 ==========
fprintf('正在保存结果...\n');

% 5.1 保存为 MAT 文件
matFile = sprintf('%sresult.mat', outputPrefix);
save(matFile, 'fc', 'bandSPL_dB', 'bandRMS', 'FreqBands', ...
     'Fs', 'nChannels', 'audioFile', 'p_ref', 'audioData');
fprintf('  MAT文件: %s\n', matFile);

% 5.2 保存为 CSV 文件
csvFile = sprintf('%sresult.csv', outputPrefix);
fid = fopen(csvFile, 'w');

% 写入表头
fprintf(fid, '中心频率(Hz),下限(Hz),上限(Hz)');
for ch = 1:nChannels
    fprintf(fid, ',通道%d_SPL(dB)', ch);
end
fprintf(fid, '\n');

% 写入数据
for i = 1:length(fc)
    fprintf(fid, '%.1f,%.2f,%.2f', fc(i), FreqBands(i,1), FreqBands(i,2));
    for ch = 1:nChannels
        fprintf(fid, ',%.4f', bandSPL_dB(i, ch));
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('  CSV文件: %s\n', csvFile);

%% ========== 6. 绘制图形 ==========
fprintf('正在绘制图形...\n');

% 6.1 主频谱图
if nChannels <= 4
    % 少通道：单图显示
    figure('Position', [100, 100, 1200, 600]);
    colors = lines(nChannels);
    hold on;
    for ch = 1:nChannels
        semilogx(fc, bandSPL_dB(:, ch), '-o', ...
                 'LineWidth', 2, 'MarkerSize', 6, ...
                 'Color', colors(ch,:), ...
                 'DisplayName', sprintf('通道 %d', ch));
    end
    grid on;
    xlabel('中心频率 (Hz)', 'FontSize', 12);
    ylabel('声压级 (dB SPL)', 'FontSize', 12);
    title(['1/3倍频程谱 - ' audioFile], 'FontSize', 14, 'Interpreter', 'none');
    legend('Location', 'best');
    xlim([fc(1)*0.9, fc(end)*1.1]);
else
    % 多通道：子图显示
    numRows = ceil(nChannels / 2);
    figure('Position', [100, 100, 1400, 300*numRows]);
    for ch = 1:nChannels
        subplot(numRows, 2, ch);
        semilogx(fc, bandSPL_dB(:, ch), '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
        grid on;
        xlabel('中心频率 (Hz)');
        ylabel('声压级 (dB SPL)');
        title(sprintf('通道 %d', ch));
        xlim([fc(1)*0.9, fc(end)*1.1]);
    end
    sgtitle(['1/3倍频程谱 - ' audioFile], 'Interpreter', 'none');
end

% 保存图形
figFile = sprintf('%splot.png', outputPrefix);
saveas(gcf, figFile);
fprintf('  图形文件: %s\n', figFile);

% 6.2 热力图（如果通道数>1）
if nChannels > 1
    figure('Position', [150, 150, 1200, 600]);
    imagesc(1:nChannels, 1:length(fc), bandSPL_dB);
    colorbar;
    colormap('jet');

    set(gca, 'YDir', 'normal');
    set(gca, 'YTick', 1:length(fc));
    set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%.1f', x), fc, 'UniformOutput', false));

    xlabel('通道号', 'FontSize', 12);
    ylabel('中心频率 (Hz)', 'FontSize', 12);
    title(['1/3倍频程谱热力图 (dB SPL) - ' audioFile], 'FontSize', 14, 'Interpreter', 'none');

    % 保存热力图
    heatmapFile = sprintf('%sheatmap.png', outputPrefix);
    saveas(gcf, heatmapFile);
    fprintf('  热力图: %s\n', heatmapFile);
end

fprintf('\n========== 计算完成 ==========\n');
fprintf('所有结果已保存到当前目录\n');
