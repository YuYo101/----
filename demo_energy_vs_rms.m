%% 演示：能量、功率、RMS的关系
% 说明为什么从能量到RMS声压需要除以N
clear; clc;

fprintf('========== 能量、功率、RMS 的关系演示 ==========\n\n');

%% 1. 生成测试信号
Fs = 1000;  % 1 kHz 采样
T = 1;      % 1 秒
t = (0:1/Fs:T-1/Fs)';
N = length(t);

% 正弦波信号：幅度 A = 2 Pa，频率 100 Hz
A = 2;
f = 100;
signal = A * sin(2*pi*f*t);

fprintf('测试信号参数:\n');
fprintf('  幅度 A = %.2f Pa\n', A);
fprintf('  频率 f = %d Hz\n', f);
fprintf('  样本数 N = %d\n\n', N);

%% 2. 理论计算
fprintf('【理论值】\n');

% 正弦波的 RMS 值
theoretical_RMS = A / sqrt(2);
fprintf('正弦波 RMS = A/√2 = %.4f Pa\n', theoretical_RMS);

% 正弦波的平均功率
theoretical_Power = (A^2) / 2;
fprintf('平均功率 = A²/2 = %.4f Pa²\n', theoretical_Power);

% 正弦波的总能量（1秒）
theoretical_Energy = theoretical_Power * N;  % 功率 × 样本数
fprintf('总能量 = 功率 × N = %.4f Pa² · samples\n\n', theoretical_Energy);

%% 3. 直接从时域计算
fprintf('【时域计算】\n');

% 能量
time_Energy = sum(signal.^2);
fprintf('时域能量 = Σx² = %.4f Pa² · samples\n', time_Energy);

% 平均功率
time_Power = time_Energy / N;
fprintf('平均功率 = 能量/N = %.4f Pa²\n', time_Power);

% RMS
time_RMS = sqrt(time_Power);
fprintf('RMS = √(功率) = %.4f Pa\n', time_RMS);

% 验证：直接计算RMS
time_RMS_direct = sqrt(mean(signal.^2));
fprintf('RMS (直接) = √(mean(x²)) = %.4f Pa\n\n', time_RMS_direct);

%% 4. 通过 computeBandEnergy 计算
fprintf('【通过 computeBandEnergy 计算】\n');

% 定义覆盖信号频率的频带
FreqBands = [50, 150];  % 50-150 Hz（包含100Hz）

% 计算频带能量
bandEnergy = computeBandEnergy(signal, Fs, FreqBands);
fprintf('computeBandEnergy 结果 = %.4f Pa² · samples\n', bandEnergy);

% 转换为功率
bandPower = bandEnergy / N;
fprintf('频带功率 = 能量/N = %.4f Pa²\n', bandPower);

% 转换为RMS
bandRMS = sqrt(bandPower);
fprintf('频带RMS = √(功率) = %.4f Pa\n', bandRMS);

% 或者一步到位
bandRMS_direct = sqrt(bandEnergy / N);
fprintf('频带RMS (一步) = √(能量/N) = %.4f Pa\n\n', bandRMS_direct);

%% 5. 对比结果
fprintf('【结果对比】\n');
fprintf('方法                  | RMS (Pa) | 误差\n');
fprintf('----------------------|----------|----------\n');
fprintf('理论值 (A/√2)         | %.4f   | -\n', theoretical_RMS);
fprintf('时域计算              | %.4f   | %.2f%%\n', time_RMS, abs(time_RMS - theoretical_RMS)/theoretical_RMS*100);
fprintf('computeBandEnergy     | %.4f   | %.2f%%\n\n', bandRMS, abs(bandRMS - theoretical_RMS)/theoretical_RMS*100);

%% 6. 关键公式总结
fprintf('========== 关键公式总结 ==========\n');
fprintf('1. 能量 E = Σx² = %.4f Pa² · samples\n', time_Energy);
fprintf('2. 功率 P = E/N = %.4f Pa²\n', time_Power);
fprintf('3. RMS = √P = √(E/N) = %.4f Pa\n', time_RMS);
fprintf('4. 声压级 SPL = 20×log₁₀(RMS/p_ref)\n\n');

fprintf('为什么要除以N？\n');
fprintf('  因为 RMS 的定义就是"均方根"：\n');
fprintf('  - 均（Mean）：除以样本数 N\n');
fprintf('  - 方（Square）：信号的平方\n');
fprintf('  - 根（Root）：开平方根\n');
fprintf('  所以 RMS = √(mean(x²)) = √(Σx²/N) = √(E/N)\n\n');

%% 7. 错误示范
fprintf('========== 如果不除以N会怎样？ ==========\n');

wrong_RMS = sqrt(bandEnergy);  % 错误：直接对能量开方
fprintf('❌ 错误做法: RMS = √(能量) = %.4f Pa\n', wrong_RMS);
fprintf('   这个值远大于真实RMS (%.4f Pa)\n', theoretical_RMS);
fprintf('   误差: %.1f倍\n\n', wrong_RMS / theoretical_RMS);

correct_RMS = sqrt(bandEnergy / N);  % 正确：先除以N
fprintf('✓ 正确做法: RMS = √(能量/N) = %.4f Pa\n', correct_RMS);
fprintf('   与理论值一致！\n\n');

fprintf('========================================\n');
