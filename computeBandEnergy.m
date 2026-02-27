function bandPower = computeBandEnergy(signal, Fs, FreqBands)
% computeBandEnergy - 计算信号在各自定义频带的平均功率
%
% 输入:
%   signal - 输入信号矩阵 [样本数 x 通道数]
%   Fs - 采样频率 (Hz)
%   FreqBands - 频带定义矩阵 [频带数 x 2], 每行为[下限频率, 上限频率]
%
% 输出:
%   bandPower - 各频带平均功率 [频带数 x 通道数] (单位: Pa²)
%
% 注意:
%   - 输出为平均功率，单位为 Pa²
%   - 如需RMS声压，使用: RMS = sqrt(bandPower)
%   - 如需总能量，使用: Energy = bandPower * N (N为样本数)
%
% 修改说明:
%   - 2026-02-27: 改为输出平均功率而非总能量（除以样本数N）
%   - 使用单边频谱避免重复计算
%   - 正确处理FFT对称性
%   - 添加Nyquist频率边界检查

[N, nChannels] = size(signal);
numBands = size(FreqBands, 1);

% 检查频带是否超出Nyquist频率
if any(FreqBands(:,2) > Fs/2)
    warning('computeBandEnergy:FreqBandExceedsNyquist', ...
            '频带上限超出Nyquist频率 (%.2f Hz)，将被截断', Fs/2);
end

bandPower = zeros(numBands, nChannels);

% FFT参数
NFFT = 2^nextpow2(N);
% 只计算单边频谱（0到Nyquist频率）
freqVector = (0:NFFT/2) * Fs / NFFT;

for ch = 1:nChannels
    % 计算FFT
    signalFFT = fft(signal(:, ch), NFFT);

    % 单边功率谱（只取前半部分，0到Nyquist）
    singleSidedFFT = signalFFT(1:NFFT/2+1);
    % Parseval定理：使用NFFT归一化，确保时域和频域能量一致
    powerSpectrum = abs(singleSidedFFT).^2 / NFFT;

    % 补偿对称性：除DC(1)和Nyquist(end)外，其他频率能量乘以2
    powerSpectrum(2:end-1) = 2 * powerSpectrum(2:end-1);

    % 计算各频带平均功率（除以样本数N）
    for band = 1:numBands
        fLow = FreqBands(band, 1);
        fHigh = min(FreqBands(band, 2), Fs/2);  % 不超过Nyquist频率

        % 找到频带范围内的频率索引（包含上边界）
        bandIdx = (freqVector >= fLow) & (freqVector <= fHigh);

        % 累加频带内的功率谱密度并除以N得到平均功率
        bandPower(band, ch) = sum(powerSpectrum(bandIdx)) / N;
    end
end

end
