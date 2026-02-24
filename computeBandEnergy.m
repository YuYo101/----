function bandEnergy = computeBandEnergy(signal, Fs, FreqBands)
% computeBandEnergy - 计算信号在各自定义频带的能量
%
% 输入:
%   signal - 输入信号矩阵 [样本数 x 通道数]
%   Fs - 采样频率 (Hz)
%   FreqBands - 频带定义矩阵 [频带数 x 2], 每行为[下限频率, 上限频率]
%
% 输出:
%   bandEnergy - 各频带能量 [频带数 x 通道数]

[N, nChannels] = size(signal);
numBands = size(FreqBands, 1);
bandEnergy = zeros(numBands, nChannels);

% FFT参数
NFFT = 2^nextpow2(N);
freqVector = (0:NFFT-1) * Fs / NFFT;

for ch = 1:nChannels
    % 计算FFT
    signalFFT = fft(signal(:, ch), NFFT);
    powerSpectrum = abs(signalFFT).^2 / N;

    % 计算各频带能量
    for band = 1:numBands
        fLow = FreqBands(band, 1);
        fHigh = FreqBands(band, 2);

        % 找到频带范围内的频率索引
        bandIdx = (freqVector >= fLow) & (freqVector < fHigh);

        % 累加频带内的能量
        bandEnergy(band, ch) = sum(powerSpectrum(bandIdx));
    end
end

end
