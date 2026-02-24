%基于gammatonefiltbank和multibandEQ的声场复现仿真
%% Step1 传递函数导入 
clc;clear;
NchS=30;NchR=10;
[Tempfile, FilePath] = uigetfile({'*.mat','mat文件 (*.mat)';...
    '*.*','所有文件 (*.*)'},'请选择灵敏度文件');
% FilePath='D:\声品质数据库\SFR\SFR_Auto\';
load([FilePath,'Sensi.mat']);
Sensi=Sensi(1:NchR);
%% Step2 目标声读取
% wav
[OriginTargetSound, Fs] = audioread([FilePath,'10ch_10s_hp30_lp20000_20260209_road_shuinilu_speed_40_03.wav']);
% OriginTargetSound=data;
TargetSound=OriginTargetSound(:,1:NchR);
N=size(OriginTargetSound,1);
T=N/Fs;
%% Step3 传函读取
tic;
load([FilePath,'ir.mat']);
ir=ir(:,1:NchR,1:NchS);
% ir=resample(ir,Fs,44100);
% 优化: 使用向量化FFT一次性计算所有通道
HestAll=zeros(N,NchR,NchS);
for i=1:NchR
    HestAll(:,i,:)=fft(ir(:,i,:), N, 1);  % 沿第一维度FFT,向量化处理
end
% 计算svd结果
% 优化: 预分配并直接使用三维矩阵S索引
U_svd=zeros(N,NchR,NchR);
S_svd=zeros(N,NchR,NchS);
V_svd=zeros(N,NchS,NchS);
for nf=1:N
    % 优化: 直接使用squeeze提取二维切片,避免双重循环
    G = squeeze(HestAll(nf,:,:));
    [U_svd(nf,:,:),S_svd(nf,:,:),V_svd(nf,:,:)] = svd(G);
end
clear HestAll G ir
endtime3=toc;
%% Step4 正则化参数
tic;
lmd0=1;
lmd = zeros(N, min(NchS,NchR));
%混合正则化参数
s0=1/2/lmd0;
% 优化: 向量化提取对角线元素和正则化参数计算
L_s=zeros(N,NchS,NchR);
L_Re=zeros(N,NchR,NchR);
for nf=1:N
    s = diag(squeeze(S_svd(nf,:,:)));
    % 向量化计算正则化参数
    lmd_temp = zeros(min(NchS,NchR),1);
    lmd_temp(s>=1/s0) = 0;
    mask_mid = (s>1/2/s0) & (s<1/s0);
    lmd_temp(mask_mid) = sqrt(s(mask_mid)/s0-s(mask_mid).^2);
    lmd_temp(s<=1/2/s0) = 1/2/s0;
    lmd(nf,:) = lmd_temp';

    % 向量化计算L_s矩阵
    S_diag = diag(squeeze(S_svd(nf,:,:)));
    for k=1:min(NchS,NchR)
        L_s(nf,k,k) = S_diag(k)/(S_diag(k)^2+lmd(nf,k)^2);
        L_Re(nf,k,k)=S_diag(k)^2/(S_diag(k)^2+lmd(nf,k)^2);
    end
end
endtime4=toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step5 驱动信号计算
tic;
TargetSound_V = TargetSound.*Sensi;  % 将声压数据转换为电压数据[Pa]*[V/Pa] = [V]
yro=TargetSound_V(1:N,1:NchR);
Ro = fft(yro,N);          % 使用fft函数将yro矩阵变为频域信号R
% 进行声场复现计算
% 优化: 预计算VLU矩阵,避免反馈循环中重复计算
S = zeros(NchS, N);
Re = zeros(NchR, N);
VLU=zeros(N,NchS,NchR);
VLU_Re=zeros(N,NchR,NchR);
for nf = 1:N                                % 采样点设置
    % 预计算VLU = V * L_s * U'
    VLU(nf,:,:)=(squeeze(V_svd(nf,:,:))) * (squeeze(L_s(nf,:,:))) * (squeeze(U_svd(nf,:,:)))';
    VLU_Re(nf,:,:)=(squeeze(U_svd(nf,:,:)))*(squeeze(L_Re(nf,:,:)))*(squeeze(U_svd(nf,:,:)))';
    S(:,nf)=(squeeze(VLU(nf,:,:)))*(Ro(nf,1:NchR).');
    Re(:,nf)=(squeeze(VLU_Re(nf,:,:)))*(Ro(nf,1:NchR).');%与S*G计算得到的结果是一致的%%%%仿真代码%%%%%
end
eest = ifft(S.');
audiowrite([FilePath, 'Reproduction.wav'],eest,Fs,"BitsPerSample",64);
Re_T=ifft(Re.');
audiowrite([FilePath, 'ReproductionResponse.wav'],Re_T,Fs,"BitsPerSample",64);
endtime5=toc;
%% Step6 误差计算
tic;
ReproductionSound=Re_T(1:N,1:NchR);
%导出Feedback变量到工作区，重置反馈声为目标声
yo_fb = TargetSound_V;   %%%循环外
%基于自定义频带的误差计算
% 定义自定义频带: 30个频带，20Hz开始，每10Hz一个频带
numcf = 30;  % 频带数量
FreqBands = zeros(numcf, 2);
for i = 1:numcf
    FreqBands(i,1) = 20 + (i-1)*10;  % 下限频率
    FreqBands(i,2) = 20 + i*10;      % 上限频率
end
FreqCenter = mean(FreqBands, 2);  % 中心频率

% 计算目标声和复现声在各频带的能量
O = computeBandEnergy(TargetSound_V(:,1:NchR), Fs, FreqBands);
R = computeBandEnergy(ReproductionSound(:,1:NchR), Fs, FreqBands);
Error = 10.*log10(O./R);
Error = Error.';
EdB=[FreqCenter.'; Error];
save([FilePath,'Error.mat'],"EdB");
FreL=25;FreU=315;             %%%%%%%反馈频段%%%%%%%%
NumL=find(FreqCenter==FreL);NumU=find(FreqCenter>=FreU,1,'first');
FBDone=0;
endtime6=toc;
%% Step7 反馈循环
tic;
Nfb=3;
for Nfbi=1+FBDone:Nfb+FBDone
    FBDone_in=0;
    %增益计算
    Gain = zeros(NchR,numcf);
    Gm = zeros(1,NchR);
    IGm = zeros(1,NchR);

    for i = NumL : NumU 
        for j = 1:NchR        
            Gain(j,i) = Error(j,i);
        end
        [Gm(1,i),IGm(1,i)] = max(abs(Gain(:,i)));
    end

    % 调整不同频率区间的各传声器增益情况
    for i = 1:length(Gain(1,:))
        [Gm(1,i),IGm(1,i)] = max(abs(Gain(:,i))); % 求出某个频率区间的最大值
        for j = 1:NchR
            if Gm(1,i)>2 && j~=IGm(1,i)
                Gain(j,i) = Gain(j,i) +...
                    sign(Gain(IGm(1,i),i))*0.2;
            end
        end
    end
    %反馈声源信号计算
    % 使用自定义频带均衡器
    for i = 1:NchR
        yo_fb(:,i) = applyCustomEQ(yo_fb(:,i), Gain(i,:), FreqBands, Fs);
    end
    % 获取多通道传声器记录的目标声场数据存放于yo矩阵中
    yo = yo_fb;     %均衡后的TargetSound_V
    Ro = fft(yo,N);           % 使用fft函数将yro矩阵变为频域信号Ro
    % 优化: 直接使用预计算的VLU矩阵,避免重复SVD和矩阵乘法
    S = zeros(NchS, N);
    FB = zeros(NchR, N);
    for nf = 1:N                            % 采样点设置
        S(:,nf)=squeeze(VLU(nf,:,:))*(Ro(nf,1:NchR).');
        FB(:,nf)=squeeze(VLU_Re(nf,:,:))*(Ro(nf,1:NchR).');
    end
    eest = ifft(S.');
    audiowrite([FilePath, ['Feedback',num2str(Nfbi),'.wav']],eest,Fs,"BitsPerSample",64);
    eest_FB = ifft(FB.');
    audiowrite([FilePath, ['FeedbackResponse',num2str(Nfbi),'.wav']],eest_FB,Fs,"BitsPerSample",64);
    FeedbackSound=eest_FB(1:N,1:NchR);
    %基于自定义频带的误差计算
    O = computeBandEnergy(TargetSound_V(:,1:NchR), Fs, FreqBands);
    R = computeBandEnergy(FeedbackSound(:,1:NchR), Fs, FreqBands);
    Error = 10.*log10(O./R);
    Error = Error.';
    EdB=[FreqCenter.'; Error];
    save([FilePath,'Error',num2str(Nfbi),'.mat'],"EdB");
    FBDone_in=FBDone_in+1;
end
FBDone=FBDone+FBDone_in;
endtime7=toc;
%% 后处理
%% PSD
nfft = round(Fs/1);
win = hann(nfft);
nov = floor(nfft *0.5);
[psd_Tar, f] = pwelch(TargetSound, win, nov, nfft, Fs);
[psd_Re, ~] = pwelch(ReproductionSound./Sensi, win, nov, nfft, Fs);
[psd_FB, ~] = pwelch(FeedbackSound./Sensi, win, nov, nfft, Fs);
psd_Tar_dB=10*log10(psd_Tar./4e-10);
psd_Re_dB=10*log10(psd_Re./4e-10);
psd_FB_dB=10*log10(psd_FB./4e-10);
% 重采样
% FsO=48000;
% nfftO = round(FsO/1);
% winO = hann(nfftO);
% novO = floor(nfftO *0.5);
% [psd_Re, fO] = pwelch(ReproductionSound, winO, novO, nfftO, FsO);
% psd_Re_dB=10*log10(psd_Re./4e-10);
for i=1:1
    % figure
    plot(f,psd_Tar_dB(:,i),'LineWidth',1,'LineStyle','-')
    hold on
    plot(f,psd_Re_dB(:,i),'LineWidth',1,'color','r','LineStyle','-.')
    hold on
    plot(f,psd_FB_dB(:,i),'LineWidth',1,'color','b','LineStyle',':')
    set(gca, 'XScale', 'line','XLim',[20,20000],'XGrid','on')
    % xticks([20,40,60,80,100,200,300,400,500,1000,2000,3000,4000,5000]);
    % legend('i8','x11')
    xlabel('频率/Hz'),ylabel('声压级/dB')
end