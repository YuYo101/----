# 自定义频带修改说明

## 修改日期
2026-02-24

## 修改内容

### 从1/3倍频程改为自定义线性频带

**原方案**:
- 使用MATLAB的poctave函数计算1/3倍频程能量
- 使用graphicEQ均衡器（固定30个中心频率: 25Hz-20000Hz）
- 频带间隔为对数分布

**新方案**:
- 自定义30个线性频带，每个频带宽度10Hz
- 频带范围: 20-30Hz, 30-40Hz, 40-50Hz, ..., 300-310Hz
- 基于FFT的频域能量计算和均衡

---

## 修改细节

### 1. 频带定义 (SFR_Auto.m:185-191)

```matlab
% 定义自定义频带: 30个频带，20Hz开始，每10Hz一个频带
numcf = 30;  % 频带数量
FreqBands = zeros(numcf, 2);
for i = 1:numcf
    FreqBands(i,1) = 20 + (i-1)*10;  % 下限频率
    FreqBands(i,2) = 20 + i*10;      % 上限频率
end
FreqCenter = mean(FreqBands, 2);  % 中心频率
```

**频带列表**:
| 频带编号 | 频率范围(Hz) | 中心频率(Hz) |
|---------|-------------|------------|
| 1 | 20-30 | 25 |
| 2 | 30-40 | 35 |
| 3 | 40-50 | 45 |
| ... | ... | ... |
| 29 | 300-310 | 305 |
| 30 | 310-320 | 315 |

### 2. 误差计算 (SFR_Auto.m:194-196)

**修改前**:
```matlab
O = poctave(TargetSound_V(...), Fs, 'FrequencyLimits', [20,20000], 'BandsPerOctave', 3);
R = poctave(ReproductionSound(...), Fs, 'FrequencyLimits', [20,20000], 'BandsPerOctave', 3);
```

**修改后**:
```matlab
O = computeBandEnergy(TargetSound_V(...), Fs, FreqBands);
R = computeBandEnergy(ReproductionSound(...), Fs, FreqBands);
```

### 3. 均衡器 (SFR_Auto.m:231-233)

**修改前**:
```matlab
equalizer = graphicEQ('EQOrder',12,'Bandwidth','1/3 octave', 'Structure','Cascade', 'SampleRate',Fs);
for i = 1:NchR
    equalizer.Gains = Gain(i,numcf-29:end);
    yo_fb(:,i) = equalizer(yo_fb(:,i));
end
```

**修改后**:
```matlab
for i = 1:NchR
    yo_fb(:,i) = applyCustomEQ(yo_fb(:,i), Gain(i,:), FreqBands, Fs);
end
```

### 4. 反馈频段调整 (SFR_Auto.m:200)

**修改前**: `FreL=20; FreU=5000;` (20Hz-5000Hz)
**修改后**: `FreL=20; FreU=300;` (20Hz-300Hz)

---

## 新增辅助函数

### computeBandEnergy.m

**功能**: 计算信号在各自定义频带的能量

**输入**:
- `signal`: 输入信号矩阵 [样本数 x 通道数]
- `Fs`: 采样频率 (Hz)
- `FreqBands`: 频带定义矩阵 [频带数 x 2]

**输出**:
- `bandEnergy`: 各频带能量 [频带数 x 通道数]

**算法**:
1. 对每个通道信号进行FFT
2. 计算功率谱密度
3. 对每个频带，累加该频带范围内的功率谱能量

### applyCustomEQ.m

**功能**: 对信号应用自定义频带均衡

**输入**:
- `inputSignal`: 输入信号 [样本数 x 1]
- `gains`: 各频带增益(dB) [频带数 x 1]
- `FreqBands`: 频带定义矩阵 [频带数 x 2]
- `Fs`: 采样频率 (Hz)

**输出**:
- `outputSignal`: 均衡后的信号 [样本数 x 1]

**算法**:
1. 对输入信号进行FFT
2. 在频域构建增益曲线（将dB转换为线性增益）
3. 对正负频率对称应用增益
4. 逆FFT回到时域

**优势**:
- 避免级联滤波器的相位失真
- 精确控制各频带增益
- 计算效率高

---

## 技术优势

### 与1/3倍频程方案对比

| 特性 | 1/3倍频程 | 自定义线性频带 |
|------|----------|---------------|
| 频带分布 | 对数分布 | 线性分布 |
| 低频分辨率 | 较低 | 高（10Hz精度） |
| 高频分辨率 | 较高 | 仅覆盖20-320Hz |
| 适用场景 | 全频段分析 | 低频精细控制 |
| 计算复杂度 | 中等 | 低 |
| 相位失真 | 有（级联滤波器） | 无（频域处理） |

### 适用场景

自定义线性频带特别适合:
- **低频噪声控制**: 车辆路噪、发动机怠速噪声
- **结构振动**: 20-300Hz的结构共振模式
- **精细调谐**: 需要精确控制特定频率的场景

---

## 使用注意事项

### 1. 频带范围限制
- 当前实现仅覆盖20-320Hz
- 如需更宽频率范围，修改频带定义:
  ```matlab
  numcf = 50;  % 增加频带数量
  bandWidth = 20;  % 改为20Hz带宽
  ```

### 2. 采样频率要求
- 确保 `Fs >= 2 * 最高频率`
- 推荐 `Fs >= 1000Hz` (当前最高频率320Hz)

### 3. 性能考虑
- FFT计算复杂度: O(N log N)
- 对于长信号，考虑分段处理

---

## Git版本管理

### 提交历史
```
* 1e46f93 (HEAD -> master) 修改为自定义频带计算(30个频带，20Hz起，每10Hz一个频带)
* 06fae8b 优化SFR_Auto.m以减少计算量
* 678da01 重置为原始版本
```

### 回退到1/3倍频程版本
```bash
git checkout 06fae8b -- SFR_Auto.m
```

### 查看修改对比
```bash
git diff 06fae8b 1e46f93
```

---

## 测试建议

1. **功能测试**: 验证误差计算和均衡效果
2. **频率响应**: 检查各频带的实际频率响应
3. **相位特性**: 确认无相位失真
4. **边界条件**: 测试极端增益值（如±20dB）

---

## 未来扩展

### 可能的改进方向

1. **自适应频带**:
   ```matlab
   % 根据信号特性自动调整频带宽度
   FreqBands = adaptiveBandDefinition(signal, Fs);
   ```

2. **平滑过渡**:
   ```matlab
   % 在频带边界使用汉宁窗平滑
   gainCurve = smoothGainTransition(gainCurve, FreqBands);
   ```

3. **非均匀频带**:
   ```matlab
   % 低频密集，高频稀疏
   FreqBands = [20 25; 25 30; 30 40; 40 60; 60 100; ...];
   ```

---

## 维护者
修改实施: Claude Code
日期: 2026-02-24
文件: SFR_Auto.m, computeBandEnergy.m, applyCustomEQ.m
