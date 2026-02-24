%%SFR_Auto
%% 参数设置
clc;clear;
NchS=30;NchR=10;  
[Tempfile, FilePath] = uigetfile({'*.mat','mat文件 (*.mat)';...
    '*.*','所有文件 (*.*)'},'请选择灵敏度文件');
% FilePath='D:\声品质数据库\SFR\SFR_Auto\';
load([FilePath,'Sensi.mat']);
Sensi=Sensi(1:NchR);
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
%% 目标声读取
% wav
% [OriginTargetSound, Fs] = audioread([FilePath,'TargetSound.wav']);
OriginTargetSound=data;
TargetSound=OriginTargetSound(:,1:NchR);
N=size(OriginTargetSound,1);
T=N/Fs;
%% 传函读取
load([FilePath,'ir.mat']);
ir=ir(:,1:NchR,1:NchS);
% ir=resample(ir,Fs,44100);
% 优化: 使用向量化FFT一次性计算所有通道
HestAll=zeros(N,NchR,NchS);
for i=1:NchR
    HestAll(:,i,:)=fft(ir(:,i,:), N, 1);  % 沿第一维度FFT,向量化处理
end
% 计算svd结果
% 优化: 预分配并直接使用三维矩阵索引
U_svd=zeros(N,NchR,NchR);
S_svd=zeros(N,NchR,NchR);
V_svd=zeros(N,NchR,NchR);
for nf=1:N
    % 优化: 直接使用squeeze提取二维切片,避免双重循环
    G = squeeze(HestAll(nf,:,:));
    [U_svd(nf,:,:),S_svd(nf,:,:),V_svd(nf,:,:)] = svd(G);
end
%% 正则化参数
lmd0=1;
lmd = zeros(N, min(NchS,NchR));
%混合正则化参数
s0=1/2/lmd0;
% 优化: 向量化提取对角线元素和正则化参数计算
L_s=zeros(N,NchS,NchR);
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
    end
end
%% 驱动信号计算
TargetSound_V = TargetSound.*Sensi;  % 将声压数据转换为电压数据[Pa]*[V/Pa] = [V]
yro=TargetSound_V(1:N,1:NchR);
Ro = fft(yro,N);          % 使用fft函数将yro矩阵变为频域信号R
% 进行声场复现计算
% 优化: 预计算VLU矩阵,避免反馈循环中重复计算
S = zeros(NchS, N);
VLU=zeros(N,NchS,NchR);
for nf = 1:N                                % 采样点设置
    % 预计算VLU = V * L_s * U'
    VLU(nf,:,:)=(squeeze(V_svd(nf,:,:))) * (squeeze(L_s(nf,:,:))) * (squeeze(U_svd(nf,:,:)))';
    S(:,nf)=(squeeze(VLU(nf,:,:)))*(Ro(nf,1:NchR).');
end
Eest = S.';
eest = ifft(Eest);
audiowrite([FilePath, 'Reproduction.wav'],eest,Fs,"BitsPerSample",64);
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
[ReproductionSound, ~] = audioread(AudioFromDevice);
ReproductionSound=ReproductionSound(1:N,1:NchR);
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
O = computeBandEnergy(TargetSound_V(BufferSize*10+1:end-BufferSize*10,1:NchR), Fs, FreqBands);
R = computeBandEnergy(ReproductionSound(BufferSize*10+1:end-BufferSize*10,1:NchR), Fs, FreqBands);
Error = 10.*log10(O./R);
Error = Error.';
EdB=[FreqCenter.'; Error];
save([FilePath,'Error.mat'],"EdB");
FreL=20;FreU=300;             %%%%%%%反馈频段%%%%%%%%
NumL=find(FreqCenter==FreL);NumU=find(FreqCenter>=FreU,1,'first');
FBDone=0;
%%%%%%%%%   反馈循环     %%%%%%%%%
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
    for nf = 1:N                            % 采样点设置
        S(:,nf)=squeeze(VLU(nf,:,:))*(Ro(nf,1:NchR).');
    end
    Eest = S.';
    eest = ifft(Eest);
    audiowrite([FilePath, ['Feedback',num2str(Nfbi),'.wav']],eest,Fs,"BitsPerSample",64);%%保存复现声源信号
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
    O = computeBandEnergy(TargetSound_V(BufferSize*10+1:end-BufferSize*10,1:NchR), Fs, FreqBands);
    R = computeBandEnergy(FeedbackSound(BufferSize*10+1:end-BufferSize*10,1:NchR), Fs, FreqBands);
    Error = 10.*log10(O./R);
    Error = Error.';
    EdB=[FreqCenter.'; Error];
    save([FilePath,'Error',num2str(Nfbi),'.mat'],"EdB");
    FBDone_in=FBDone_in+1;
end
FBDone=FBDone+FBDone_in;
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