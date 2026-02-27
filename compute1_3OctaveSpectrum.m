function [fc, bandSPL_dB, bandRMS, FreqBands] = compute1_3OctaveSpectrum(audioData, Fs, options)
% compute1_3OctaveSpectrum - 计算音频信号的1/3倍频程谱（声压级 dB SPL）
%
% 语法:
%   [fc, bandSPL_dB, bandRMS, FreqBands] = compute1_3OctaveSpectrum(audioData, Fs)
%   [fc, bandSPL_dB, bandRMS, FreqBands] = compute1_3OctaveSpectrum(audioData, Fs, options)
%
% 输入:
%   audioData - 音频信号矩阵 [样本数 x 通道数]
%   Fs        - 采样频率 (Hz)
%   options   - 可选参数结构体（Name-Value对）:
%               'FreqRange'  - 频率范围 [fmin, fmax] (默认: [20, 20000] Hz)
%               'RefPressure'- 参考声压 (默认: 2e-5 Pa = 20 μPa)
%
% 输出:
%   fc         - 中心频率向量 [频带数 x 1] (Hz)
%   bandSPL_dB - 各频带声压级 [频带数 x 通道数] (dB SPL)
%   bandRMS    - 各频带RMS声压 [频带数 x 通道数] (Pa)
%   FreqBands  - 频带边界 [频带数 x 2]，每行为 [下限, 上限] (Hz)
%
% 算法说明:
%   1. 使用标准1/3倍频程中心频率（IEC 61260标准）
%   2. 通过 computeBandEnergy 计算各频带平均功率
%   3. 转换为RMS声压: RMS = sqrt(平均功率)
%   4. 计算声压级: SPL = 20 * log10(RMS / p_ref)
%
% 示例:
%   % 读取音频文件并计算1/3倍频程谱
%   [audioData, Fs] = audioread('test.wav');
%   [fc, SPL] = compute1_3OctaveSpectrum(audioData, Fs);
%
%   % 绘制频谱
%   semilogx(fc, SPL);
%   xlabel('频率 (Hz)'); ylabel('声压级 (dB SPL)');
%   grid on;
%
% 参考:
%   IEC 61260:1995 - Electroacoustics - Octave-band and fractional-octave-band filters
%
% 作者: Claude Code
% 日期: 2026-02-27
% 修改: 2026-02-27 - 修正声压级计算公式，使用正确的参考声压

%% 参数解析
arguments
    audioData (:,:) double
    Fs (1,1) double {mustBePositive}
    options.FreqRange (1,2) double {mustBePositive} = [20, 20000]
    options.RefPressure (1,1) double {mustBePositive} = 2e-5
end

[N, nChannels] = size(audioData);
p_ref = options.RefPressure;  % 参考声压 (Pa)

%% 定义1/3倍频程频带
% 标准1/3倍频程中心频率（IEC 61260标准）
fc_standard = [20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, ...
               250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, ...
               2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000];

% 筛选有效频率范围
fc_valid_mask = (fc_standard >= options.FreqRange(1)) & ...
                (fc_standard <= options.FreqRange(2)) & ...
                (fc_standard <= Fs/2);  % 不超过Nyquist频率
fc = fc_standard(fc_valid_mask)';

% 计算频带边界
% 1/3倍频程带宽比 = 2^(1/6) ≈ 1.122
bandwidth_ratio = 2^(1/6);
fl = fc / bandwidth_ratio;  % 下限频率
fu = fc * bandwidth_ratio;  % 上限频率

FreqBands = [fl, fu];
numBands = length(fc);

%% 计算各频带功率
% 使用 computeBandEnergy 函数（现在返回平均功率）
bandPower = computeBandEnergy(audioData, Fs, FreqBands);

%% 转换为声压级
% 计算RMS声压（均方根声压）
% RMS = sqrt(平均功率)
bandRMS = sqrt(bandPower);  % [频带数 x 通道数]

% 计算声压级（dB SPL）
% SPL = 20 * log10(p_rms / p_ref)
bandSPL_dB = 20 * log10(bandRMS / p_ref + eps);  % eps防止log(0)

%% 输出提示信息
if nargout == 0
    % 如果没有输出参数，打印结果并绘图
    fprintf('1/3倍频程谱计算完成\n');
    fprintf('  采样频率: %d Hz\n', Fs);
    fprintf('  信号长度: %d 样本 (%.2f 秒)\n', N, N/Fs);
    fprintf('  通道数: %d\n', nChannels);
    fprintf('  频带数: %d\n', numBands);
    fprintf('  频率范围: %.1f - %.1f Hz\n', fc(1), fc(end));
    fprintf('  参考声压: %.2e Pa (20 μPa)\n', p_ref);
    fprintf('  SPL范围: %.2f - %.2f dB SPL\n', min(bandSPL_dB(:)), max(bandSPL_dB(:)));

    % 绘制频谱
    figure('Position', [100, 100, 1000, 600]);
    if nChannels <= 4
        % 少于4个通道：单图显示
        colors = lines(nChannels);
        hold on;
        for ch = 1:nChannels
            semilogx(fc, bandSPL_dB(:, ch), '-o', ...
                     'LineWidth', 1.5, 'MarkerSize', 5, ...
                     'Color', colors(ch,:), ...
                     'DisplayName', sprintf('通道 %d', ch));
        end
        grid on;
        xlabel('中心频率 (Hz)', 'FontSize', 12);
        ylabel('声压级 (dB SPL)', 'FontSize', 12);
        title('1/3倍频程谱', 'FontSize', 14);
        legend('Location', 'best');
    else
        % 多通道：子图显示
        numRows = ceil(nChannels / 2);
        for ch = 1:nChannels
            subplot(numRows, 2, ch);
            semilogx(fc, bandSPL_dB(:, ch), '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
            grid on;
            xlabel('中心频率 (Hz)');
            ylabel('声压级 (dB SPL)');
            title(sprintf('通道 %d', ch));
        end
    end
end

end
