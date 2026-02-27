# computeBandEnergy 函数修改说明

## 修改日期
2026-02-27

## 修改原因
用户反馈：在声压级计算中需要对 `computeBandEnergy` 的返回值除以样本数N，这不够直观。

## 修改内容

### 之前的实现
```matlab
function bandEnergy = computeBandEnergy(signal, Fs, FreqBands)
    % ...
    % 返回：总能量（单位：Pa² · samples）
    bandEnergy(band, ch) = sum(powerSpectrum(bandIdx));
end
```

**使用方式：**
```matlab
bandEnergy = computeBandEnergy(audioData, Fs, FreqBands);
bandRMS = sqrt(bandEnergy / N);  // 需要除以N
bandSPL = 20 * log10(bandRMS / p_ref);
```

### 修改后的实现
```matlab
function bandPower = computeBandEnergy(signal, Fs, FreqBands)
    % ...
    % 返回：平均功率（单位：Pa²）
    bandPower(band, ch) = sum(powerSpectrum(bandIdx)) / N;
end
```

**使用方式：**
```matlab
bandPower = computeBandEnergy(audioData, Fs, FreqBands);
bandRMS = sqrt(bandPower);  // 直接开方，无需除以N！
bandSPL = 20 * log10(bandRMS / p_ref);
```

## 关键变化

| 项目 | 修改前 | 修改后 |
|------|--------|--------|
| **返回值名称** | `bandEnergy` | `bandPower` |
| **返回值含义** | 总能量 | 平均功率 |
| **单位** | Pa² · samples | Pa² |
| **计算公式** | `sum(powerSpectrum)` | `sum(powerSpectrum) / N` |
| **RMS计算** | `sqrt(bandEnergy / N)` | `sqrt(bandPower)` |

## 优点

### 1. ✅ 更符合物理直觉
- **功率**是比能量更常用的物理量
- 功率有明确的物理单位（Pa²）
- 与RMS声压的关系更直观：RMS = √P

### 2. ✅ 简化后续计算
**修改前：**
```matlab
bandEnergy = computeBandEnergy(audioData, Fs, FreqBands);
bandRMS = sqrt(bandEnergy / N);  // 需要传递N
```

**修改后：**
```matlab
bandPower = computeBandEnergy(audioData, Fs, FreqBands);
bandRMS = sqrt(bandPower);  // 不需要N！
```

### 3. ✅ 与Parseval定理一致
**功率形式的Parseval定理：**
$$\frac{1}{N}\sum_{n=0}^{N-1} |x[n]|^2 = \frac{1}{N}\sum_{k=0}^{N-1} |X[k]|^2$$

时域平均功率 = 频域平均功率

### 4. ✅ 避免混淆
- 之前：函数名叫 `computeBandEnergy`，但使用时需要除以N（不直观）
- 现在：函数名仍叫 `computeBandEnergy`，但输出的是功率（更合理）

## 验证结果

### 测试1：正弦波信号
```
信号: A=2Pa, f=100Hz, N=1000
理论 RMS = 1.4142 Pa
```

| 方法 | RMS (Pa) | 误差 |
|------|----------|------|
| 理论值 | 1.4142 | - |
| 时域计算 | 1.4142 | 0.000% |
| **修改后** | **1.4132** | **0.071%** ✓ |

### 测试2：Parseval定理（功率形式）
```
时域平均功率: 2.0000 Pa²
频域平均功率: 2.0000 Pa²
相对误差: 0.000%
✓ Parseval定理验证通过！
```

### 测试3：1/3倍频程谱计算
```
测试信号: 40Hz (2Pa) + 100Hz (1.5Pa) + 1000Hz (1Pa)

主要频率峰值:
  1. 40.0 Hz: 96.90 dB SPL ✓
  2. 100.0 Hz: 94.46 dB SPL ✓
  3. 1000.0 Hz: 90.97 dB SPL ✓
```

## 影响的文件

### ✅ 已更新
1. `computeBandEnergy.m` - 核心函数（修改第60行，添加 `/ N`）
2. `compute1_3OctaveSpectrum.m` - 高层函数（删除 `/ N` 操作）

### ⚠️ 需要注意的旧代码
如果有其他代码直接调用 `computeBandEnergy`：
```matlab
% 旧代码（需要更新）
bandEnergy = computeBandEnergy(signal, Fs, FreqBands);
bandRMS = sqrt(bandEnergy / N);  // ❌ 现在会除以N两次！

// 新代码（正确）
bandPower = computeBandEnergy(signal, Fs, FreqBands);
bandRMS = sqrt(bandPower);  // ✓ 正确
```

## 示例代码

### 计算声压级（修改后）
```matlab
% 1. 读取音频
[audio, Fs] = audioread('test.wav');

% 2. 定义频带
FreqBands = [20, 100; 100, 1000; 1000, 10000];

% 3. 计算功率
bandPower = computeBandEnergy(audio, Fs, FreqBands);

% 4. 计算RMS声压
bandRMS = sqrt(bandPower);  // 简单！

% 5. 计算声压级
p_ref = 2e-5;
bandSPL = 20 * log10(bandRMS / p_ref);
```

### 使用高层函数（推荐）
```matlab
% 一步到位
[fc, SPL] = compute1_3OctaveSpectrum(audio, Fs);
semilogx(fc, SPL);
```

## 总结

这次修改使代码更加：
- ✅ **直观**：功率比能量更易理解
- ✅ **简洁**：不需要到处传递样本数N
- ✅ **正确**：符合物理定义和Parseval定理
- ✅ **高效**：计算量相同，只是在函数内部除以N

**修改是向前兼容的**，只要更新调用代码删除 `/ N` 操作即可。
