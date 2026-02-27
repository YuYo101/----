# 频带功率分析方法对比：FFT vs Welch

## 概述

本代码库提供两种频带功率计算方法：
1. **computeBandEnergy** - 基于单次FFT的快速方法
2. **computeBandPower_pwelch** - 基于Welch方法的高精度方法

两者接口完全一致，可以根据应用场景灵活选择。

---

## 方法对比

### computeBandEnergy (FFT方法)

**原理**：
- 单次FFT变换
- 矩形窗（无窗）
- 直接求和频域功率

**优点**：
- ✅ **速度快**：比pwelch快8-10倍（0.025秒 vs 0.24秒）
- ✅ **实时性好**：适合在线处理
- ✅ **适合宽频带**：对于1/3倍频程等宽频带，精度足够
- ✅ **简单高效**：代码简洁，易于理解

**缺点**：
- ❌ **频率泄漏大**：矩形窗导致主瓣宽、旁瓣高
- ❌ **噪声方差大**：单次估计，方差未降低
- ❌ **噪声底高**：在噪声区域比pwelch高10-30 dB

**算法**：
```matlab
% 1. 单次FFT
signalFFT = fft(signal, NFFT);

% 2. 计算单边功率谱
singleSidedFFT = signalFFT(1:NFFT/2+1);
powerSpectrum = abs(singleSidedFFT).^2 / NFFT;
powerSpectrum(2:end-1) = 2 * powerSpectrum(2:end-1);  % DC和Nyquist除外

% 3. 频带求和
bandPower(band) = sum(powerSpectrum(bandIdx)) / N;
```

---

### computeBandPower_pwelch (Welch方法)

**原理**：
- 信号分段
- Hanning窗函数
- 分段FFT后平均（Welch方法）
- 积分功率谱密度（PSD）

**优点**：
- ✅ **频率泄漏小**：Hanning窗主瓣窄、旁瓣低
- ✅ **噪声方差小**：多段平均，大幅降低方差
- ✅ **精度高**：在信号频率处误差 < 0.1 dB
- ✅ **噪声底低**：噪声区域比FFT低10-30 dB
- ✅ **适合弱信号**：能更好地检测淹没在噪声中的信号

**缺点**：
- ❌ **计算稍慢**：比FFT慢8-10倍
- ❌ **频率分辨率损失**：窗函数展宽主瓣

**算法**：
```matlab
% 1. 使用pwelch计算PSD
[Pxx, f_psd] = pwelch(signal, window, noverlap, nfft, Fs);

% 2. 积分PSD得到功率
df = f_psd(2) - f_psd(1);  % 频率分辨率
for band = 1:numBands
    bandIdx = (f_psd >= fLow) & (f_psd <= fHigh);
    bandPower(band) = sum(Pxx(bandIdx)) * df;  % 矩形积分
end
```

**可调参数**：
```matlab
bandPower = computeBandPower_pwelch(signal, Fs, FreqBands, ...
    'Window', hanning(8192), ...  % 窗函数或窗长度
    'Overlap', 4096, ...          % 重叠样本数
    'NFFT', 16384);               % FFT点数
```

---

## 实测对比结果

### 测试1：多频信号（40, 100, 1000 Hz）

| 中心频率 | FFT方法 | pwelch方法 | 差异 |
|---------|--------|-----------|-----|
| **40 Hz** | 96.93 dB | 96.99 dB | **-0.06 dB** ✓ |
| **100 Hz** | 94.47 dB | 94.49 dB | **-0.02 dB** ✓ |
| **1000 Hz** | 90.97 dB | 90.97 dB | **-0.00 dB** ✓ |
| 20 Hz (噪声) | 55.65 dB | 38.89 dB | **+16.76 dB** |
| 63 Hz (噪声) | 54.31 dB | 43.22 dB | **+11.10 dB** |

**结论**：在信号主频率处，两者结果高度一致（< 0.1 dB）；在噪声区域，pwelch噪声底更低。

### 测试2：1Hz带宽窄带分析（100 Hz纯音）

| 频带 | FFT方法 | pwelch方法 | 差异 |
|-----|--------|-----------|-----|
| 98-99 Hz | 0.00508 Pa² | 0.00001 Pa² | +26.05 dB |
| **99-100 Hz** | **1.12678 Pa²** | **1.03318 Pa²** | **+0.38 dB** ✓ |
| **100-101 Hz** | **0.85278 Pa²** | **0.96626 Pa²** | **-0.54 dB** ✓ |
| 101-102 Hz | 0.00486 Pa² | 0.00001 Pa² | +25.59 dB |

**结论**：在主频率频带（99-100、100-101 Hz），两者误差 < 1 dB；在远离主频率的频带，FFT泄漏严重。

### 测试3：计算速度（10秒音频，1/3倍频程）

| 方法 | 耗时 | 相对速度 |
|-----|-----|---------|
| computeBandEnergy | 0.025秒 | 基准 |
| computeBandPower_pwelch | 0.239秒 | **9.6倍慢** |

**结论**：FFT方法快约10倍，适合实时处理。

---

## 使用建议

### 何时使用 computeBandEnergy（FFT方法）

✅ **推荐场景**：
- 实时处理（速度要求高）
- 宽频带分析（如1/3倍频程、1/1倍频程）
- 信噪比高的信号（主频率能量远大于噪声）
- 快速预览和原型开发
- 多通道同步处理

❌ **不适合场景**：
- 弱信号检测（信号淹没在噪声中）
- 需要精确噪声底测量
- 窄带分析（< 1Hz带宽）
- 科研级高精度测量

**典型应用**：
- 车内噪声在线监测
- 语音信号实时分析
- 音质快速评估
- 生产线噪声筛查

---

### 何时使用 computeBandPower_pwelch（Welch方法）

✅ **推荐场景**：
- 需要高精度测量（< 0.5 dB误差）
- 噪声分析（低噪声底要求）
- 窄带分析（1Hz或更窄）
- 弱信号检测
- 科研级分析和报告
- 后处理场景（非实时）

❌ **不适合场景**：
- 实时处理（速度瓶颈）
- 极短信号（< 1秒，分段不足）
- 快速预览

**典型应用**：
- 建筑声学精密测量
- 机械故障诊断（轴承振动等）
- 环境噪声监测报告
- 产品噪声认证测试
- 科研论文数据

---

## 参数调整建议

### computeBandPower_pwelch 参数调优

#### 1. Window（窗函数/窗长度）

**默认值**：`hanning(min(4096, N))`

**调整建议**：
```matlab
% 提高低频精度（增大窗长度）
bandPower = computeBandPower_pwelch(signal, Fs, FreqBands, ...
    'Window', hanning(8192));

% 处理时变信号（减小窗长度）
bandPower = computeBandPower_pwelch(signal, Fs, FreqBands, ...
    'Window', hanning(2048));

% 使用其他窗函数
bandPower = computeBandPower_pwelch(signal, Fs, FreqBands, ...
    'Window', hamming(4096));  % Hamming窗
```

**经验法则**：
- 窗长度 = Fs / (最低频率 / 2)  # 至少覆盖2个周期
- 低频分析（< 50 Hz）：窗长度 ≥ 8192
- 中频分析（50-1000 Hz）：窗长度 = 4096（默认）
- 高频分析（> 1000 Hz）：窗长度 = 2048

#### 2. Overlap（重叠样本数）

**默认值**：50%窗长度

**调整建议**：
```matlab
% 提高平滑度（增大重叠）
bandPower = computeBandPower_pwelch(signal, Fs, FreqBands, ...
    'Window', hanning(4096), 'Overlap', 3072);  % 75%重叠

% 加快速度（减小重叠）
bandPower = computeBandPower_pwelch(signal, Fs, FreqBands, ...
    'Window', hanning(4096), 'Overlap', 1024);  % 25%重叠
```

**经验法则**：
- 标准分析：50%重叠（默认）
- 高精度/平滑谱：75%重叠
- 快速计算：25%重叠
- 不推荐：< 25%（方差过大）或 > 75%（收益递减）

#### 3. NFFT（FFT点数）

**默认值**：`max(4096, 2^nextpow2(N))`

**调整建议**：
```matlab
% 提高频率分辨率（增大NFFT）
bandPower = computeBandPower_pwelch(signal, Fs, FreqBands, ...
    'NFFT', 16384);

% 宽带分析（使用默认值）
bandPower = computeBandPower_pwelch(signal, Fs, FreqBands);
```

**经验法则**：
- NFFT应该 ≥ 窗长度（否则会零填充）
- 频率分辨率 Δf = Fs / NFFT
- 精细谱分析：NFFT = 16384 或更大
- 宽带分析：NFFT = 4096（默认）

---

## 代码示例

### 基本使用（两种方法）

```matlab
% 读取音频
[signal, Fs] = audioread('test.wav');

% 定义1/3倍频程频带
fc_octave = [63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000];
bandwidth_ratio = 2^(1/6);
fl = fc_octave / bandwidth_ratio;
fu = fc_octave * bandwidth_ratio;
FreqBands = [fl', fu'];

%% 方法1：快速FFT方法
bandPower_FFT = computeBandEnergy(signal, Fs, FreqBands);

%% 方法2：高精度pwelch方法
bandPower_pwelch = computeBandPower_pwelch(signal, Fs, FreqBands);

%% 计算声压级
bandRMS_FFT = sqrt(bandPower_FFT);
bandSPL_FFT = 20 * log10(bandRMS_FFT / 2e-5);

bandRMS_pwelch = sqrt(bandPower_pwelch);
bandSPL_pwelch = 20 * log10(bandRMS_pwelch / 2e-5);

%% 绘制对比
figure;
semilogx(fc_octave, bandSPL_FFT, 'b-o', 'DisplayName', 'FFT方法');
hold on;
semilogx(fc_octave, bandSPL_pwelch, 'r--s', 'DisplayName', 'pwelch方法');
xlabel('中心频率 (Hz)');
ylabel('声压级 (dB SPL)');
title('1/3倍频程谱对比');
legend;
grid on;
```

### 高级使用（自定义参数）

```matlab
% 场景：低频精密测量（20-200 Hz）
fc_low = [20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200];
fl_low = fc_low / 2^(1/6);
fu_low = fc_low * 2^(1/6);
FreqBands_low = [fl_low', fu_low'];

% 使用长窗、高重叠、高FFT分辨率
bandPower_lowfreq = computeBandPower_pwelch(signal, Fs, FreqBands_low, ...
    'Window', hanning(8192), ...   % 长窗提高低频精度
    'Overlap', 6144, ...           % 75%重叠提高平滑度
    'NFFT', 32768);                % 高FFT分辨率

% 转换为声压级
bandSPL_lowfreq = 20 * log10(sqrt(bandPower_lowfreq) / 2e-5);
```

### 多通道处理

```matlab
% 读取多通道音频（如10通道）
[audioData, Fs] = audioread('10ch_audio.wav');  % [N × 10]

% 定义频带
fc_octave = [20, 25, 31.5, ..., 20000];
% ... (定义FreqBands)

% 一次性计算所有通道
bandPower_all = computeBandPower_pwelch(audioData, Fs, FreqBands);
% bandPower_all: [numBands × 10]

% 提取各通道结果
bandSPL_ch1 = 20 * log10(sqrt(bandPower_all(:, 1)) / 2e-5);
bandSPL_ch2 = 20 * log10(sqrt(bandPower_all(:, 2)) / 2e-5);
% ...

% 绘制多通道对比
figure;
for ch = 1:10
    SPL = 20 * log10(sqrt(bandPower_all(:, ch)) / 2e-5);
    semilogx(fc_octave, SPL, 'DisplayName', sprintf('通道%d', ch));
    hold on;
end
xlabel('频率 (Hz)');
ylabel('声压级 (dB SPL)');
title('10通道1/3倍频程谱');
legend;
grid on;
```

---

## 完整文件列表

### 核心函数
1. **computeBandEnergy.m** - FFT方法（快速）
2. **computeBandPower_pwelch.m** - pwelch方法（高精度）
3. **compute1_3OctaveSpectrum.m** - 高层函数（1/3倍频程）

### 对比脚本
4. **compare_two_methods.m** - 主要对比脚本（推荐先运行）
5. **compare_1Hz_corrected.m** - 1Hz带宽对比（验证窄带一致性）
6. **compare_computeBandEnergy_vs_pwelch.m** - 详细对比

### 使用示例
7. **example_use_pwelch_method.m** - pwelch方法使用指南（6个示例）
8. **quick_octave_analysis.m** - 快速分析脚本
9. **example_use_octave_spectrum.m** - 1/3倍频程使用示例

### 验证脚本
10. **verify_computeBandEnergy.m** - Parseval定理验证
11. **demo_energy_vs_rms.m** - 能量/功率/RMS概念说明

### 文档
12. **CHANGELOG_computeBandEnergy.md** - 修改日志
13. **README_octave_analysis.md** - 1/3倍频程分析文档
14. **README_frequency_band_analysis.md** - 本文档

---

## 常见问题（FAQ）

### Q1: 两种方法能得到相同结果吗？

**A**: 取决于场景：
- **信号主频率处**（信噪比高）：差异 < 0.1 dB，可认为一致 ✓
- **噪声区域**（信噪比低）：差异 10-30 dB，pwelch更准确
- **1Hz窄带**：在主频率差异 < 1 dB，远离主频率差异很大

**建议**：如果只关心主要频率成分，两种方法均可；如果关心噪声底，必须用pwelch。

### Q2: 为什么pwelch在信号频率处反而比FFT略小？

**A**: 主要原因：
1. **窗函数损失**：Hanning窗会损失一定能量（约1.36 dB）
2. **主瓣展宽**：窗函数展宽主瓣，能量分散到相邻频率
3. **分段截断**：信号被分段，边界效应

但这种差异极小（< 0.5 dB），在工程精度范围内。

### Q3: 如何选择窗长度？

**A**: 基本原则：
- **低频分析**：窗长度 = Fs / (最低频率 / 2)，至少2个周期
  - 例如：20 Hz → 窗长度 ≥ 44100/10 = 4410 样本
- **时变信号**：窗长度不宜过长，以捕捉时变特性
- **折中方案**：4096是常见默认值（Fs=44.1kHz时约0.09秒）

### Q4: 哪种方法更"正确"？

**A**: 没有绝对正确，取决于应用：
- **FFT方法**：更接近"瞬时能量"，适合宽带快速评估
- **pwelch方法**：更接近"统计平均"，适合精密测量

声学测量标准（如IEC 61672）通常基于FFT的思路，但科研分析更青睐pwelch。

### Q5: 为什么函数叫computeBandEnergy但输出是功率？

**A**: 历史原因。该函数最初返回"总能量"（单位：Pa²·samples），需要除以样本数N才能得到功率。后来为了简化使用，修改为直接输出"平均功率"（单位：Pa²），但函数名保留了。详见 `CHANGELOG_computeBandEnergy.md`。

**正确理解**：
```matlab
bandPower = computeBandEnergy(...);  % 虽然叫Energy，实际返回Power
bandRMS = sqrt(bandPower);           % 直接开方即可，无需再除以N
```

---

## 理论背景

### Parseval定理（功率形式）

$$\frac{1}{N}\sum_{n=0}^{N-1} |x[n]|^2 = \frac{1}{N}\sum_{k=0}^{N-1} |X[k]|^2$$

时域平均功率 = 频域平均功率

### 单边功率谱

对于实信号，FFT具有共轭对称性，因此：

$$P_{\text{single-sided}}(k) = \begin{cases}
|X(0)|^2 / N & k = 0 \\
2|X(k)|^2 / N & 0 < k < N/2 \\
|X(N/2)|^2 / N & k = N/2
\end{cases}$$

### Welch方法

信号分段，每段加窗、FFT、计算PSD，最后平均：

$$\hat{P}_{xx}(f) = \frac{1}{K} \sum_{i=1}^{K} \left| \text{FFT}(x_i \cdot w) \right|^2$$

方差降低因子 ≈ K（段数），但频率分辨率下降。

### 频带功率积分

频带功率 = PSD在频带内的积分：

$$P_{\text{band}} = \int_{f_L}^{f_H} S_{xx}(f) \, df \approx \sum_{f=f_L}^{f_H} S_{xx}(f) \cdot \Delta f$$

---

## 性能基准测试

**测试环境**：
- CPU: (取决于运行环境)
- MATLAB: R2021a或更高
- 信号：10秒，44.1 kHz，单通道
- 频带：21个1/3倍频程（20-2000 Hz）

**结果**：

| 方法 | 时间 | 相对速度 | 内存 |
|-----|------|---------|------|
| computeBandEnergy | 0.025秒 | 1× | ~5 MB |
| computeBandPower_pwelch (默认) | 0.239秒 | 9.6× | ~15 MB |
| computeBandPower_pwelch (8192窗) | 0.450秒 | 18× | ~20 MB |

**结论**：FFT方法快10倍，pwelch方法内存占用多2-3倍。

---

## 引用与参考

1. Welch, P. (1967). "The use of fast Fourier transform for the estimation of power spectra: A method based on time averaging over short, modified periodograms". *IEEE Transactions on Audio and Electroacoustics*, 15(2), 70-73.

2. IEC 61260:1995. "Electroacoustics - Octave-band and fractional-octave-band filters".

3. MATLAB Documentation: `pwelch`, `fft`, `bandpower`.

4. Oppenheim, A. V., & Schafer, R. W. (2009). *Discrete-Time Signal Processing* (3rd ed.). Prentice Hall.

---

## 更新日志

- **2026-02-27**: 创建`computeBandPower_pwelch`函数，完成两种方法对比测试
- **2026-02-27**: 修改`computeBandEnergy`从输出能量改为输出功率
- **2026-02-27**: 创建完整的对比脚本和使用示例

---

## 联系与反馈

如有问题或建议，请通过以下方式联系：
- 代码仓库：[待添加]
- 邮箱：[待添加]

---

**最后更新**：2026-02-27
**作者**：Claude Code
**版本**：1.0
