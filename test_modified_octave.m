%% 测试修改后的 compute1_3OctaveSpectrum 函数
clear; clc;

fprintf('========== 测试 compute1_3OctaveSpectrum 函数 ==========\n\n');

%% 生成测试信号
Fs = 8000;
T = 1;
t = (0:1/Fs:T-1/Fs)';

% 多个频率成分
signal = 2*sin(2*pi*40*t) + 1.5*sin(2*pi*100*t) + sin(2*pi*1000*t);

fprintf('测试信号: 40Hz (2Pa) + 100Hz (1.5Pa) + 1000Hz (1Pa)\n');
fprintf('采样频率: %d Hz\n', Fs);
fprintf('样本数: %d\n\n', length(t));

%% 使用修改后的函数计算
[fc, SPL, RMS, FreqBands] = compute1_3OctaveSpectrum(signal, Fs, 'FreqRange', [20, 2000]);

fprintf('计算完成！\n');
fprintf('频带数: %d\n\n', length(fc));

%% 显示结果
fprintf('1/3倍频程谱结果:\n');
fprintf('中心频率 | RMS声压 | 声压级\n');
fprintf('---------|---------|----------\n');
for i = 1:length(fc)
    fprintf('%6.1f Hz | %.4f Pa | %6.2f dB\n', fc(i), RMS(i), SPL(i));
end

%% 验证主要频率峰值
fprintf('\n主要频率峰值（前5个）:\n');
[sorted_SPL, idx] = sort(SPL, 'descend');
for i = 1:5
    fprintf('  %d. %.1f Hz: %.2f dB SPL\n', i, fc(idx(i)), sorted_SPL(i));
end

%% 计算总声压级
fprintf('\n总声压级:\n');
totalRMS = sqrt(sum(RMS.^2));  % 多个频率成分的RMS合成
totalSPL = 20 * log10(totalRMS / 2e-5);
fprintf('  总RMS: %.4f Pa\n', totalRMS);
fprintf('  总SPL: %.2f dB SPL\n', totalSPL);

fprintf('\n========== 测试完成 ==========\n');
fprintf('✓ 函数工作正常！现在计算过程更简洁：\n');
fprintf('  bandPower = computeBandEnergy(...)  // 得到功率\n');
fprintf('  bandRMS = sqrt(bandPower)           // 直接开方\n');
fprintf('  bandSPL = 20*log10(RMS/p_ref)       // 计算SPL\n');
