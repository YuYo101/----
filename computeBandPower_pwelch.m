function bandPower = computeBandPower_pwelch(signal, Fs, FreqBands, options)
% computeBandPower_pwelch - 基于pwelch方法计算信号在各频带的平均功率
%
% 语法:
%   bandPower = computeBandPower_pwelch(signal, Fs, FreqBands)
%   bandPower = computeBandPower_pwelch(signal, Fs, FreqBands, options)
%
% 输入:
%   signal    - 输入信号矩阵 [样本数 x 通道数]
%   Fs        - 采样频率 (Hz)
%   FreqBands - 频带定义矩阵 [频带数 x 2], 每行为[下限频率, 上限频率]
%   options   - 可选参数结构体（Name-Value对）:
%               'Window'   - 窗函数或窗长度 (默认: hanning(min(4096, N)))
%               'Overlap'  - 重叠样本数 (默认: 50%窗长度)
%               'NFFT'     - FFT点数 (默认: max(4096, 2^nextpow2(N)))
%
% 输出:
%   bandPower - 各频带平均功率 [频带数 x 通道数] (单位: Pa²)
%
% 与 computeBandEnergy 的区别:
%   computeBandEnergy:
%     - 使用单次FFT + 矩形窗
%     - 快速，适合实时计算
%     - 频率泄漏较大
%   computeBandPower_pwelch:
%     - 使用pwelch方法（分段平均 + Hanning窗）
%     - 精度更高，频率泄漏更小
%     - 噪声方差更小
%     - 计算稍慢
%
% 算法说明:
%   1. 使用pwelch计算功率谱密度（PSD）
%   2. 对每个频带，积分PSD: Power = ∫ PSD(f) df
%   3. 积分采用矩形法: Power ≈ Σ(PSD × Δf)
%
% 示例:
%   % 基本使用
%   [signal, Fs] = audioread('test.wav');
%   FreqBands = [20, 100; 100, 1000; 1000, 10000];
%   bandPower = computeBandPower_pwelch(signal, Fs, FreqBands);
%   bandRMS = sqrt(bandPower);
%   bandSPL = 20 * log10(bandRMS / 2e-5);
%
%   % 自定义窗函数
%   bandPower = computeBandPower_pwelch(signal, Fs, FreqBands, ...
%                                       'Window', hanning(8192));
%
% 参考:
%   - Welch, P. (1967). "The use of fast Fourier transform for the
%     estimation of power spectra"
%   - MATLAB pwelch 文档
%
% 作者: Claude Code
% 日期: 2026-02-27

%% 参数解析
arguments
    signal (:,:) double
    Fs (1,1) double {mustBePositive}
    FreqBands (:,2) double
    options.Window = []
    options.Overlap = []
    options.NFFT = []
end

[N, nChannels] = size(signal);
numBands = size(FreqBands, 1);

% 检查频带是否超出Nyquist频率
if any(FreqBands(:,2) > Fs/2)
    warning('computeBandPower_pwelch:FreqBandExceedsNyquist', ...
            '频带上限超出Nyquist频率 (%.2f Hz)，将被截断', Fs/2);
end

%% 设置默认参数
% 窗函数
if isempty(options.Window)
    window_length = min(4096, N);
    window = hanning(window_length);
else
    if isscalar(options.Window)
        % 如果是标量，表示窗长度
        window = hanning(options.Window);
    else
        % 否则直接使用提供的窗
        window = options.Window(:);
    end
end
window_length = length(window);

% 重叠
if isempty(options.Overlap)
    noverlap = floor(window_length / 2);  % 默认50%重叠
else
    noverlap = options.Overlap;
end

% FFT点数
if isempty(options.NFFT)
    nfft = max(4096, 2^nextpow2(N));
else
    nfft = options.NFFT;
end

%% 初始化输出
bandPower = zeros(numBands, nChannels);

%% 对每个通道计算
for ch = 1:nChannels
    % 使用pwelch计算功率谱密度
    [Pxx, f_psd] = pwelch(signal(:, ch), window, noverlap, nfft, Fs);

    % 频率分辨率
    df = f_psd(2) - f_psd(1);

    % 对每个频带积分PSD
    for band = 1:numBands
        fLow = FreqBands(band, 1);
        fHigh = min(FreqBands(band, 2), Fs/2);  % 不超过Nyquist频率

        % 找到频带范围内的频率索引（包含上边界）
        bandIdx = (f_psd >= fLow) & (f_psd <= fHigh);

        % 积分PSD得到功率（矩形积分法）
        % Power = ∫ PSD(f) df ≈ Σ(PSD × Δf)
        bandPower(band, ch) = sum(Pxx(bandIdx)) * df;
    end
end

end
