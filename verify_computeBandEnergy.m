%% 验证 computeBandEnergy.m 函数的正确性
% 测试内容：
% 1. Parseval定理验证：频域能量 = 时域能量
% 2. 多频率成分信号测试
% 3. 边界条件测试

clear; clc;
fprintf('========== computeBandEnergy 验证测试 ==========\n\n');

%% 测试1: Parseval定理验证
fprintf('【测试1】Parseval定理验证\n');
fprintf('----------------------------------------------\n');

% 生成测试信号
Fs = 1000;  % 采样频率 1000 Hz
T = 1;      % 信号时长 1 秒
t = (0:1/Fs:T-1/Fs)';
N = length(t);

% 多频率成分信号（3个通道）
f1 = 10;   % 10 Hz
f2 = 50;   % 50 Hz
f3 = 150;  % 150 Hz
A1 = 1; A2 = 2; A3 = 1.5;

signal = zeros(N, 3);
signal(:,1) = A1*sin(2*pi*f1*t) + A2*cos(2*pi*f2*t);
signal(:,2) = A2*sin(2*pi*f2*t) + A3*cos(2*pi*f3*t);
signal(:,3) = A1*sin(2*pi*f1*t) + A3*sin(2*pi*f3*t);

% 计算时域总能量
timeEnergy = sum(signal.^2);

% 定义覆盖全频段的频带
FreqBands = [0, Fs/2];  % 从0到Nyquist频率

% 计算频域能量
bandEnergy = computeBandEnergy(signal, Fs, FreqBands);
freqEnergy = bandEnergy;

% 比较时域和频域能量
fprintf('通道 | 时域能量 | 频域能量 | 相对误差\n');
fprintf('-----|----------|----------|----------\n');
for ch = 1:3
    relError = abs(timeEnergy(ch) - freqEnergy(ch)) / timeEnergy(ch) * 100;
    fprintf('  %d  | %8.4f | %8.4f | %6.3f%%\n', ...
            ch, timeEnergy(ch), freqEnergy(ch), relError);
end

% 判断Parseval定理是否满足（误差<1%）
if all(abs(timeEnergy - freqEnergy) ./ timeEnergy < 0.01)
    fprintf('✓ Parseval定理验证通过！\n\n');
else
    fprintf('✗ Parseval定理验证失败！\n\n');
end

%% 测试2: 多频带分解测试
fprintf('【测试2】多频带能量分解\n');
fprintf('----------------------------------------------\n');

% 定义多个频带
FreqBands = [
    0,   20;   % 低频带（包含10Hz）
    20,  100;  % 中频带（包含50Hz）
    100, 200;  % 高频带（包含150Hz）
    200, Fs/2  % 超高频带（无信号）
];

% 计算各频带能量
bandEnergy = computeBandEnergy(signal, Fs, FreqBands);

% 显示结果
fprintf('频带范围 | 通道1能量 | 通道2能量 | 通道3能量\n');
fprintf('---------|-----------|-----------|----------\n');
bandNames = {'0-20 Hz  ', '20-100 Hz', '100-200Hz', '200-500Hz'};
for band = 1:size(FreqBands, 1)
    fprintf('%s | %9.4f | %9.4f | %9.4f\n', ...
            bandNames{band}, bandEnergy(band,1), bandEnergy(band,2), bandEnergy(band,3));
end

% 验证各频带能量之和 = 总能量
totalBandEnergy = sum(bandEnergy, 1);
fprintf('\n各频带能量之和:\n');
fprintf('通道 | 总能量 | 时域能量 | 相对误差\n');
fprintf('-----|--------|----------|----------\n');
for ch = 1:3
    relError = abs(timeEnergy(ch) - totalBandEnergy(ch)) / timeEnergy(ch) * 100;
    fprintf('  %d  | %6.4f | %8.4f | %6.3f%%\n', ...
            ch, totalBandEnergy(ch), timeEnergy(ch), relError);
end

if all(abs(timeEnergy - totalBandEnergy) ./ timeEnergy < 0.01)
    fprintf('✓ 多频带分解验证通过！\n\n');
else
    fprintf('✗ 多频带分解验证失败！\n\n');
end

%% 测试3: 频率成分定位验证
fprintf('【测试3】频率成分定位验证\n');
fprintf('----------------------------------------------\n');

% 通道1应该在0-20Hz和20-100Hz有能量
% 通道2应该在20-100Hz和100-200Hz有能量
% 通道3应该在0-20Hz和100-200Hz有能量

fprintf('预期能量分布:\n');
fprintf('通道1: 主要在 0-20Hz (10Hz成分) 和 20-100Hz (50Hz成分)\n');
fprintf('通道2: 主要在 20-100Hz (50Hz成分) 和 100-200Hz (150Hz成分)\n');
fprintf('通道3: 主要在 0-20Hz (10Hz成分) 和 100-200Hz (150Hz成分)\n\n');

% 检查能量占比
for ch = 1:3
    fprintf('通道%d各频带能量占比:\n', ch);
    energyRatio = bandEnergy(:,ch) / sum(bandEnergy(:,ch)) * 100;
    for band = 1:size(FreqBands, 1)
        fprintf('  %s: %5.2f%%\n', bandNames{band}, energyRatio(band));
    end
    fprintf('\n');
end

%% 测试4: 边界条件测试
fprintf('【测试4】边界条件测试\n');
fprintf('----------------------------------------------\n');

% 测试超出Nyquist频率的频带
fprintf('测试超出Nyquist频率的频带定义...\n');
FreqBands_exceed = [0, 600];  % 超过Fs/2=500Hz

% 捕获警告
warning('off', 'computeBandEnergy:FreqBandExceedsNyquist');
lastwarn('');
bandEnergy_exceed = computeBandEnergy(signal(:,1), Fs, FreqBands_exceed);
[warnMsg, warnId] = lastwarn;

if strcmp(warnId, 'computeBandEnergy:FreqBandExceedsNyquist')
    fprintf('✓ 成功检测到超出Nyquist频率并发出警告\n');
else
    fprintf('✗ 未能检测到超出Nyquist频率\n');
end
warning('on', 'computeBandEnergy:FreqBandExceedsNyquist');

% 验证能量是否仍然正确
relError = abs(timeEnergy(1) - bandEnergy_exceed) / timeEnergy(1) * 100;
fprintf('边界截断后能量误差: %.3f%%\n', relError);

if relError < 1
    fprintf('✓ 边界截断处理正确\n\n');
else
    fprintf('✗ 边界截断处理有误\n\n');
end

%% 测试5: 单频信号理论验证
fprintf('【测试5】单频信号理论验证\n');
fprintf('----------------------------------------------\n');

% 生成纯正弦波
f_test = 100;  % 100 Hz
A_test = 2;    % 幅度为2
signal_sine = A_test * sin(2*pi*f_test*t);

% 理论能量：E = A^2 * N / 2（对于正弦波）
theory_energy = A_test^2 * N / 2;

% 定义窄频带围绕100Hz
FreqBands_narrow = [
    0,   90;
    90,  110;
    110, Fs/2
];

bandEnergy_sine = computeBandEnergy(signal_sine, Fs, FreqBands_narrow);

fprintf('频带范围 | 能量值 | 能量占比\n');
fprintf('---------|--------|----------\n');
total_energy_sine = sum(bandEnergy_sine);
for band = 1:3
    ratio = bandEnergy_sine(band) / total_energy_sine * 100;
    fprintf('%s | %6.2f | %6.2f%%\n', ...
            ['[' num2str(FreqBands_narrow(band,1)) '-' num2str(FreqBands_narrow(band,2)) 'Hz]'], ...
            bandEnergy_sine(band), ratio);
end

fprintf('\n时域能量: %.4f\n', sum(signal_sine.^2));
fprintf('频域能量: %.4f\n', total_energy_sine);
fprintf('理论能量: %.4f\n', theory_energy);

if bandEnergy_sine(2) / total_energy_sine > 0.99
    fprintf('✓ 单频信号能量集中在正确频带\n');
else
    fprintf('✗ 单频信号能量分布异常\n');
end

%% 总结
fprintf('\n========== 验证测试完成 ==========\n');
fprintf('所有测试项均已完成，请检查上述结果。\n');
fprintf('如果所有测试都显示 ✓，则函数工作正常。\n');
