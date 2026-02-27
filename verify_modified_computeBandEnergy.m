%% 验证修改后的 computeBandEnergy 函数（现在输出功率）
clear; clc;

fprintf('========== 验证修改后的 computeBandEnergy 函数 ==========\n\n');

%% 1. 生成测试信号
Fs = 1000;  % 1 kHz
T = 1;      % 1 秒
t = (0:1/Fs:T-1/Fs)';
N = length(t);

% 正弦波：幅度 A = 2 Pa，频率 100 Hz
A = 2;
f = 100;
signal = A * sin(2*pi*f*t);

fprintf('测试信号:\n');
fprintf('  幅度: %.2f Pa\n', A);
fprintf('  频率: %d Hz\n', f);
fprintf('  样本数: %d\n\n', N);

%% 2. 理论值
theoretical_RMS = A / sqrt(2);
theoretical_Power = theoretical_RMS^2;

fprintf('理论值:\n');
fprintf('  RMS声压 = A/√2 = %.4f Pa\n', theoretical_RMS);
fprintf('  平均功率 = RMS² = %.4f Pa²\n\n', theoretical_Power);

%% 3. 使用修改后的 computeBandEnergy 计算
fprintf('使用修改后的 computeBandEnergy:\n');

FreqBands = [50, 150];  % 50-150 Hz（包含100Hz）
bandPower = computeBandEnergy(signal, Fs, FreqBands);

fprintf('  返回值（平均功率）= %.4f Pa²\n', bandPower);

% 计算RMS（现在不需要除以N了！）
bandRMS = sqrt(bandPower);
fprintf('  RMS = √(功率) = %.4f Pa\n', bandRMS);

% 计算声压级
p_ref = 2e-5;
bandSPL = 20 * log10(bandRMS / p_ref);
fprintf('  声压级 = %.2f dB SPL\n\n', bandSPL);

%% 4. 与时域直接计算对比
fprintf('时域直接计算:\n');
time_RMS = sqrt(mean(signal.^2));
time_Power = time_RMS^2;
time_SPL = 20 * log10(time_RMS / p_ref);

fprintf('  RMS = %.4f Pa\n', time_RMS);
fprintf('  平均功率 = %.4f Pa²\n', time_Power);
fprintf('  声压级 = %.2f dB SPL\n\n', time_SPL);

%% 5. 验证结果
fprintf('【验证结果】\n');
fprintf('方法                   | RMS (Pa) | 功率(Pa²) | 误差\n');
fprintf('-----------------------|----------|-----------|--------\n');
fprintf('理论值 (A/√2)          | %.4f   | %.4f    | -\n', theoretical_RMS, theoretical_Power);
fprintf('时域计算               | %.4f   | %.4f    | %.3f%%\n', ...
        time_RMS, time_Power, abs(time_RMS - theoretical_RMS)/theoretical_RMS*100);
fprintf('computeBandEnergy      | %.4f   | %.4f    | %.3f%%\n', ...
        bandRMS, bandPower, abs(bandRMS - theoretical_RMS)/theoretical_RMS*100);

%% 6. 关键改进说明
fprintf('\n========== 关键改进 ==========\n');
fprintf('【修改前】\n');
fprintf('  bandEnergy = computeBandEnergy(...)  // 返回总能量\n');
fprintf('  bandRMS = sqrt(bandEnergy / N)       // 需要除以N\n\n');

fprintf('【修改后】\n');
fprintf('  bandPower = computeBandEnergy(...)   // 返回平均功率\n');
fprintf('  bandRMS = sqrt(bandPower)            // 直接开方即可！\n\n');

fprintf('优点:\n');
fprintf('  1. ✓ 更符合物理直觉（功率是常用物理量）\n');
fprintf('  2. ✓ 简化后续计算（不需要传递N）\n');
fprintf('  3. ✓ 单位明确（Pa²）\n');
fprintf('  4. ✓ 与 Parseval 定理一致（时域平均功率 = 频域平均功率）\n');

%% 7. 验证 Parseval 定理（使用功率）
fprintf('\n========== Parseval 定理验证（功率形式） ==========\n');

% 计算覆盖全频段的功率
FreqBands_full = [0, Fs/2];
freq_Power = computeBandEnergy(signal, Fs, FreqBands_full);

% 时域平均功率
time_Power_avg = mean(signal.^2);

fprintf('时域平均功率: %.4f Pa²\n', time_Power_avg);
fprintf('频域平均功率: %.4f Pa²\n', freq_Power);
fprintf('相对误差: %.3f%%\n', abs(time_Power_avg - freq_Power)/time_Power_avg*100);

if abs(time_Power_avg - freq_Power)/time_Power_avg < 0.01
    fprintf('✓ Parseval定理验证通过（功率形式）！\n');
else
    fprintf('✗ Parseval定理验证失败！\n');
end

fprintf('\n========== 验证完成 ==========\n');
