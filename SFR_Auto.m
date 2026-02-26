%%SFR_Auto
%% Step1 参数设置 
clc;clear;
NchS=30;NchR=10;
InputChannelMap=[1:NchR];
FilePath = uigetdir;
% FilePath='D:\声品质数据库\SFR\SFR_Auto\';
load([FilePath,'\Sensi.mat']);
Sensi=Sensi(InputChannelMap);
%% 目标声导入
%%连接数据库
dbHost = "nvh-dmu-test-default.rdsm4xq6y4k2pkl.rds.bj.baidubce.com";
dbPort = 13307;
dbName = "nvh_dmu";
dbUser = "tanghao_rw";
dbPass = "#BpJgVkZ5VIxjse8)Y8";

javaaddpath("D:\声品质数据库\data\env\mysql-connector-j-8.0.33.jar");

conn = database(dbName, dbUser, dbPass, ...
    "Vendor","MySQL", ...
    "Server",dbHost, ...
    "PortNumber",dbPort);
%%%
if isopen(conn)
    disp("数据库连接成功！");
    %%查询索引%%
    % sql = 'select data_id, major, problem_condition, measure_point from sound_quality_time_data_th';
    % data0 = fetch(conn, sql);
    %%读取数据%%
    sql_read = 'select * from sound_quality_time_data_th where vehicle_info="L8" and problem_condition="粗糙路-60kph" and direction = "S_in"';
    data_read = fetch(conn, sql_read);
    Fs=data_read.sample_frequency(1);
    Time=eval(data_read.x_start{1}):eval(data_read.x_step{1}):...
        (eval(data_read.x_start{1})+eval(data_read.x_step{1})*(data_read.data_length-1));
    Time=Time.';
    N=data_read.data_length(1);
    Nch=size(data_read,1);
    data=zeros(N,Nch);
    for i=1:Nch
        data(:,i)=eval(eval(data_read.y_values{i})).';
    end
else
    error("数据库连接失败");
end
% 关闭连接
close(conn);
audiowrite([FilePath,'TargetSound.wav'],data,Fs,BitsPerSample=64);
% clear conn data_read
%% Step2 目标声读取
% wav
[OriginTargetSound, Fs] = audioread([FilePath,'\10ch_10s_hp30_lp20000_20260209_road_shuinilu_speed_40_03.wav']);
% OriginTargetSound=data;
TargetSound=OriginTargetSound(:,InputChannelMap);
N=size(OriginTargetSound,1);
T=N/Fs;
%% Step3 传函读取
tic;
try
    load([FilePath,'\FRF\SVD.mat']);
    if size(S_svd,1) ~=N
        load([FilePath,'\FRF\ir.mat']);
        ir=ir(:,InputChannelMap,1:NchS);
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
        save([FilePath,'\FRF\SVD.mat'],"U_svd","S_svd","V_svd",'-v7.3','-nocompression');
        clear HestAll G ir
    end
catch
    load([FilePath,'\FRF\ir.mat']);
    ir=ir(:,InputChannelMap,1:NchS);
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
    save([FilePath,'\FRF\SVD.mat'],"U_svd","S_svd","V_svd",'-v7.3','-nocompression');
    clear HestAll G ir
end
endtime3=toc;
%% Step4 正则化参数
tic;
lmd0_base=1;  % 基准正则化参数
lmd = zeros(N, min(NchS,NchR));
% 优化: 向量化提取对角线元素和正则化参数计算
L_s=zeros(N,NchS,NchR);
s=zeros(min(NchS,NchR),N);

% 方案A: 频率依赖的正则化参数（保持共轭对称性）/能量阈值
% for nf=1:N
%     % 计算当前频点对应的实际频率 (Hz)
%     if nf <= floor(N/2)+1
%         % 正频率部分: 0 到 Fs/2
%         freq = (nf-1) * Fs / N;
%     else
%         % 负频率部分: 使用镜像点的频率（保证对称）
%         mirror_nf = N - nf + 2;
%         freq = (mirror_nf-1) * Fs / N;
%     end
% 
%     % 根据频率调整正则化参数
%     if freq < 70  % 70Hz以下，强烈放松正则化
%         lmd0 = 0.01;  % 降低10倍
%     elseif freq < 300  % 100-300Hz，中度放松
%         lmd0 = 0.3;  % 降低约3倍
%     elseif freq < 500  % 300-500Hz，轻度放松
%         lmd0 = 0.6;
%     else  % 500Hz以上，保持原值
%         lmd0 = lmd0_base;
%     end
% 
%     % 混合正则化参数计算
%     s0=1/2/lmd0;
%     s(:,nf) = diag(squeeze(S_svd(nf,:,:)));
% 
%     for i=1:min(NchS,NchR)
%         if s(i,nf)>=1/s0
%             lmd(nf,i)=0;
%         elseif s(i,nf)>1/2/s0
%             lmd(nf,i)=sqrt(s(i,nf)/s0-s(i,nf)^2);
%         else
%             lmd(nf,i)=1/2/s0;
%         end
%     end
% 
%     for k=1:min(NchS,NchR) %%奇异值矩阵
%         L_s(nf,k,k)=S_svd(nf,k,k)/(S_svd(nf,k,k)^2+lmd(nf,k)^2);
%     end
% end
%%%%%%%%%%%%%%%%
% 方案B: 固定能量传递率L_Re的正则化参数
lmd0=1;lmd_energy=0.5;
for nf=1:N
    s(:,nf) = diag(squeeze(S_svd(nf,:,:)));

    for i=1:min(NchS,NchR)
        if s(i,nf)>=lmd0
            lmd(nf,i)=0;
        else
            lmd(nf,i)=sqrt((1/lmd_energy-1)*s(i,nf)^2);
        end
    end

    for k=1:min(NchS,NchR) %%奇异值矩阵
        L_s(nf,k,k)=S_svd(nf,k,k)/(S_svd(nf,k,k)^2+lmd(nf,k)^2);
    end
end
clear S_svd
endtime4=toc;
%% Step5 驱动信号计算
tic;
TargetSound_V = TargetSound.*Sensi;  % 将声压数据转换为电压数据[Pa]*[V/Pa] = [V]
yro=TargetSound_V;
Ro = fft(yro,N);          % 使用fft函数将yro矩阵变为频域信号R
% 进行声场复现计算
% 优化: 预计算VLU矩阵,避免反馈循环中重复计算
S = zeros(NchS, N);
Re = zeros(NchR, N);
VLU=zeros(N,NchS,NchR);
for nf = 1:N                                % 采样点设置
    % 预计算VLU = V * L_s * U'
    VLU(nf,:,:)=(squeeze(V_svd(nf,:,:))) * (squeeze(L_s(nf,:,:))) * (squeeze(U_svd(nf,:,:)))';
    S(:,nf)=(squeeze(VLU(nf,:,:)))*(Ro(nf,:).');
end
eest = ifft(S.');
audiowrite([FilePath, '\Reproduction.wav'],eest,Fs,"BitsPerSample",64);
clear U_svd V_svd
endtime5=toc;
%% 播放
SafeRMS = 1;SafeLevel = 0.55;
BufferSize=2048;
% 声卡驱动
WriterDriver='ASIO';SelectedWriter='Antelope ';
ReaderDriver='ASIO';SelectedReader='';
% 选择合适的输出设备并播放复现声场信号
AudioToDevice = [FilePath, 'Reproduction.wav'];
AudioFromDevice = [FilePath, 'ReproductionResponse.wav'];
fileReader = dsp.AudioFileReader(AudioToDevice,'SamplesPerFrame',BufferSize);
fileWriter = dsp.AudioFileWriter(AudioFromDevice,'SampleRate',Fs,...
    'FileFormat','WAV','DataType','double');

deviceWriter = audioDeviceWriter(Fs,'Driver',WriterDriver,'BitDepth',...
    '24-bit integer','Device', SelectedWriter,...
    'SupportVariableSizeInput',1,'BufferSize',BufferSize);
deviceReader = audioDeviceReader(Fs,'Driver',ReaderDriver,'BitDepth',...
    '24-bit integer','Device', SelectedReader,...
    'SamplesPerFrame',BufferSize,'NumChannels',NchR);

setup(deviceWriter,zeros(BufferSize, NchS));
setup(deviceReader);
%%%%%%%%%%%  播放  %%%%%%%%%%%%%%%%%%%
% 帧数设置
totalframe = ceil(T*Fs/BufferSize); % 总帧数
framecount = 1; % 帧数计数器
% 缓入设置
NWin=10;
win = hann(2*NWin*BufferSize);
while ~isDone(fileReader)
    % 生成播放帧
    stream = fileReader();
    % 缓入
    if framecount <= NWin
        stream = win((framecount-1)*BufferSize+1:framecount*BufferSize).*stream;
    end
    % 缓出
    if framecount > totalframe - NWin
        stream = win((2*NWin-totalframe+framecount-1)*BufferSize+1: ...
            (2*NWin-totalframe+framecount)*BufferSize).*stream;
    end
    % 定义输出帧
    streamOut = stream; % 输出帧赋值，按照音频通道数调用扬声器数
    % 限制最大输出幅值及功率
    if max(max(abs(streamOut))) > SafeLevel
        break;
    end
    if max(rms(streamOut)) > SafeRMS
        break;
    end
    % 播放音频数据流
    deviceWriter(streamOut);
    audioIn = deviceReader();
    fileWriter(audioIn);
    % 计算帧数
    percent = framecount/totalframe;
    framecount = framecount + 1;
end % end while (~isDone(fileReader))
% 释放音频模块
release(fileReader);
release(deviceWriter);
release(deviceReader);
release(fileWriter);
%%%%%%%%%后处理%%%%%%%%%
%% Step6 误差计算
tic;
[ReproductionSound, ~] = audioread(AudioFromDevice);
ReproductionSound=ReproductionSound(1:N,1:NchR);
%导出Feedback变量到工作区，重置反馈声为目标声
yo_fb = TargetSound_V;   %%%循环外
%基于自定义频带的误差计算
% 定义自定义频带: 
% 20-50Hz:30个频带，20Hz开始，每1Hz一个频带
% 20-50Hz:300个频带，50Hz开始，每10Hz一个频带
numcf = 3000;  % 频带数量
FreBandwise=1; %20-50频率带宽
% FreBandwise2=10; %50-3050频率带宽
FreqBands = zeros(numcf, 2);
for i = 1:numcf
    FreqBands(i,1) = 20 + (i-1)*FreBandwise;  % 下限频率
    FreqBands(i,2) = 20 + i*FreBandwise;      % 上限频率
end
FreqCenter = mean(FreqBands, 2);  % 中心频率

% 计算目标声和复现声在各频带的能量
O = computeBandEnergy(TargetSound_V, Fs, FreqBands);
R = computeBandEnergy(ReproductionSound, Fs, FreqBands);
Error = 10.*log10(O./R);
Error = Error.';
EdB=[FreqCenter.'; Error];
save([FilePath,'\Error.mat'],"EdB");
FreL=20.5;FreU=3000;             %%%%%%%反馈频段%%%%%%%%
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
    Gm = zeros(1,NchR);%%%%%维度不对
    IGm = zeros(1,NchR);

    Gain(:,NumL:NumU)=Error(:,NumL:NumU);

    % 调整不同频率区间的各传声器增益情况
    % 方案1：非对称调整系数 - 正误差小系数，负误差增大系数
    for i = 1:length(Gain(1,:))
        [Gm(1,i),IGm(1,i)] = max(abs(Gain(:,i)));
        for j = 1:NchR
            if Gm(1,i)>2 && j~=IGm(1,i)
                if Gain(IGm(1,i),i) > 0
                    adjustStep=0.2;
                else
                    adjustStep=0.4;
                end
                Gain(j,i) = Gain(j,i) + sign(Gain(IGm(1,i),i))*adjustStep;
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
    for nf = 1:N                            % 采样点设置
        S(:,nf)=squeeze(VLU(nf,:,:))*(Ro(nf,:).');
    end
    eest = ifft(S.');
    audiowrite([FilePath, ['\Feedback',num2str(Nfbi),'.wav']],eest,Fs,"BitsPerSample",64);
    %重新播放
    AudioToDevice = [FilePath,'Feedback',num2str(Nfbi),'.wav'];
    AudioFromDevice = [FilePath,'FeedbackResponse',num2str(Nfbi),'.wav'];
    fileReader = dsp.AudioFileReader(AudioToDevice,'SamplesPerFrame',BufferSize);
    fileWriter = dsp.AudioFileWriter(AudioFromDevice,'SampleRate',Fs,...
        'FileFormat','WAV','DataType','double');

    deviceWriter = audioDeviceWriter(Fs,'Driver',WriterDriver,'BitDepth',...
        '24-bit integer','Device', SelectedWriter,...
        'SupportVariableSizeInput',1,'BufferSize',BufferSize);
    deviceReader = audioDeviceReader(Fs,'Driver',ReaderDriver,'BitDepth',...
        '24-bit integer','Device', SelectedReader,...
        'SamplesPerFrame',BufferSize,'NumChannels',NchR);
    setup(deviceWriter,zeros(BufferSize, NchS));
    setup(deviceReader);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 帧数设置
    totalframe = ceil(T*Fs/BufferSize); % 总帧数
    framecount = 1; % 帧数计数器
    % 缓入设置
    win = hann(2*NWin*BufferSize);
    while ~isDone(fileReader)
        % 取消进程
        % 生成播放帧
        stream = fileReader();
        % 缓入
        if framecount <= NWin
            stream = win((framecount-1)*BufferSize+1:framecount*BufferSize).*stream;
        end
        % 缓出
        if framecount > totalframe - NWin
            stream = win((2*NWin-totalframe+framecount-1)*BufferSize+1: ...
                (2*NWin-totalframe+framecount)*BufferSize).*stream;
        end
        % 定义输出帧
        streamOut = stream; % 输出帧赋值，按照音频通道数调用扬声器数
        % 限制最大输出幅值及功率
        if max(max(abs(streamOut))) > SafeLevel
            break;
        end
        if max(rms(streamOut)) > SafeRMS
            break;
        end
        % 播放音频数据流
        deviceWriter(streamOut);
        audioIn = deviceReader();
        fileWriter(audioIn);
        % 计算帧数
        percent = framecount/totalframe;
        framecount = framecount + 1;
    end % end while (~isDone(fileReader))
    % 释放音频模块
    release(fileReader);
    release(deviceWriter);
    release(deviceReader);
    release(fileWriter);
    [FeedbackSound, ~] = audioread(AudioFromDevice);
    FeedbackSound=FeedbackSound(1:N,1:NchR);
    %基于自定义频带的误差计算
    O = computeBandEnergy(TargetSound_V, Fs, FreqBands);
    R = computeBandEnergy(FeedbackSound, Fs, FreqBands);
    Error = 10.*log10(O./R);
    Error = Error.';
    EdB=[FreqCenter.'; Error];
    save([FilePath,'\Error',num2str(Nfbi),'.mat'],"EdB");
    FBDone_in=FBDone_in+1;
end
FBDone=FBDone+FBDone_in;
endtime7=toc;
%% 驱动信号导出至数据库
dbHost = "nvh-dmu-test-default.rdsm4xq6y4k2pkl.rds.bj.baidubce.com";
dbPort = 13307;
dbName = "nvh_dmu";
dbUser = "tanghao_rw";
dbPass = "#BpJgVkZ5VIxjse8)Y8";

javaaddpath("D:\声品质数据库\data\env\mysql-connector-j-8.0.33.jar");

conn = database(dbName, dbUser, dbPass, ...
    "Vendor","MySQL", ...
    "Server",dbHost, ...
    "PortNumber",dbPort);
%%%
if isopen(conn)
    disp("数据库连接成功！");
    %%查询索引%%
    % sql = 'select data_id, major, problem_condition, measure_point from sound_quality_time_data_th';
    % data0 = fetch(conn, sql);
    %%写入数据%%
    [audiofile,fs]=audioread('Feedback3.wav');    
    for i=1:size(audiofile,2)
    audio_json = jsonencode(audiofile(:,i));
    audio_json = jsonencode(audio_json);
    datain(i,:)=table({'路噪'},{'水泥路-40kph'},...
        {'/'},{'L8'},{'驱动信号'},...
        {['Point',num2str(i)]},{'S'},{fs}, ...
        {size(audiofile,1)},{size(audiofile,1)/fs},{'Time'},...
        {'s'},{'Pressure'},{'Pa'},...
        {audio_json},{num2str(1/fs)},{num2str(1/fs)},...
        {'于泳'},...
        'VariableNames',...
        {'major','problem_condition',...
        'problem','vehicle_info','remark',...
        'measure_point','direction','sample_frequency',...
        'data_length','time_length','x_quantity',...
        'x_unit','y_quantity','y_unit',...
        'y_values','x_start','x_step',...
        'upload_person'});
    end
    sqlwrite(conn,'sound_quality_time_data_th',datain);
else
    error("数据库连接失败");
end
% 关闭连接
close(conn);