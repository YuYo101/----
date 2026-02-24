function outputSignal = applyCustomEQ(inputSignal, gains, FreqBands, Fs)
% applyCustomEQ - 对信号应用自定义频带均衡
%
% 输入:
%   inputSignal - 输入信号 [样本数 x 1]
%   gains - 各频带增益(dB) [频带数 x 1]
%   FreqBands - 频带定义矩阵 [频带数 x 2], 每行为[下限频率, 上限频率]
%   Fs - 采样频率 (Hz)
%
% 输出:
%   outputSignal - 均衡后的信号 [样本数 x 1]

N = length(inputSignal);
numBands = size(FreqBands, 1);

% FFT参数
NFFT = 2^nextpow2(N);
freqVector = (0:NFFT-1) * Fs / NFFT;

% 计算FFT
signalFFT = fft(inputSignal, NFFT);

% 初始化频域增益曲线
gainCurve = ones(NFFT, 1);

% 对各频带应用增益
for band = 1:numBands
    fLow = FreqBands(band, 1);
    fHigh = FreqBands(band, 2);

    % 将dB转换为线性增益
    linearGain = 10^(gains(band)/20);

    % 找到频带范围内的频率索引
    bandIdx = (freqVector >= fLow) & (freqVector < fHigh);

    % 应用增益（包括正频率和负频率）
    gainCurve(bandIdx) = linearGain;

    % 对称处理负频率部分
    negBandIdx = (freqVector >= (Fs - fHigh)) & (freqVector < (Fs - fLow));
    gainCurve(negBandIdx) = linearGain;
end

% 在频域应用增益
signalFFT_EQ = signalFFT .* gainCurve;

% 逆FFT回到时域
outputSignal = real(ifft(signalFFT_EQ));

% 截取到原始长度
outputSignal = outputSignal(1:N);

end
