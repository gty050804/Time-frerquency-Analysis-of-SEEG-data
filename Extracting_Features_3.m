% 目前，运行出的gamma_epoch中包含部分Trigger（约占总数的5%）附近0.2s（由参数extent决定）的seeg数据决定，
% 第一维度为时间，第二维度为Trial，第三维度为Chn。
% 测试是否有显著响应的模块得到保留，可以更改平均值的倍数调整阈值设定。

% SessionNumber不具备鲁棒性，由于本实验每组实验次数均为2因此具备鲁棒性。

B=0;count=0;con=0;clear gamma_epoch

% 矩阵还未搭建好，完成了Trigger的分割，可以判断部分通道的响应情况。

% Selected_Chn：所选通道的响应情况

% 将Gamma_epoch的有效时间放置在每个刺激开始前的2s到开始后的2s，取决于extent

g1 = size(targetSubjects,2);   

sessionNum = 2;

Gamma_epoch_cell = cell(g1,sessionNum);   % 创建一个包含g1*2个3维矩阵的Cell，每个三维矩阵代表一个P&Session

index_of_stimulus_onset_cell = cell(length(targetSubjects),sessionNum);



%% Constrct gamma_epoch

% sel_Chn = 20;0
extent1 = 3;     % 决定Trigger前多长的时间（单位：s）
extent2 = 6;     % 决定Trigger后多长的时间（单位：s）


for subjId = targetSubjects

    % con代表了第几个P

    count = 0; con = con+1;

    fprintf("Constructing gamma-epoch of P %d \n", subjId);

    subInfo = get_subject_info(subjId);

    Fs = subInfo.Fs;   actualFs = Fs;

    sessionNo = subInfo.Session_num;

    for cons = sessionNo

        clear gamma_epoch
        
    
        fprintf("Constructing gamma-epoch of Session %d \n", cons);
    
        N = size(Subjectcell{con,cons},1);   % 当前采样点个数
    
        datacell = Subjectcell{con,cons};    % 当前预处理后的待研究矩阵
        
        Trig = datacell(:, end);
        
        index_of_stimulus_onset = find(Trig > 0);

        index_of_stimulus_onset_cell{con,cons} = index_of_stimulus_onset;
    
        % selected_Chn = datacell(:,sel_Chn);
    
        % seegDiff_smooth = smooth(abs(selected_Chn), 0.025*actualFs);
    
        g2 = size(datacell,2)-4;     % 有效通道数
    
        g3 = size(index_of_stimulus_onset,1);   % Trigger所代表的刺激数
    
        g4 = (extent1+extent2) * actualFs;               % 特定时间内的采样点数，4=2+2
    
        % gamma_epoch = zeros(g4 , g3 , g2);      % 维度1：g2->有效通道   g3->刺激试次   g4->刺激前后的一段时间
    
        for chn = 1:g2
    
            
    
            fprintf("Channel %d/%d is being processed.\n",chn,g2);
    
    
    
            selected_Chn = datacell(:,chn);
    
            seegDiff_smooth = smooth(abs(selected_Chn), 0.025*actualFs);
    
            for trial = 1:g3  % the segment number (i th)
    
                    count = count + 1;
    
                    
        
                        % start_idx = index_of_stimulus_onset(trial);
                        % 
                        % end_idx = min(start_idx + 5*actualFs, length(selected_Chn));
                        % 
                        % seegSegment = selected_Chn(start_idx:end_idx); % 5s-long EMG data segment after the trigger signal
                        % 
                        % % detect event
                        % alarm = envelop_hilbert_v2(seegSegment, round(0.025*actualFs), 1, round(0.05*actualFs), 0);
                        % 
                        % alarm_indices = find(alarm == 1);
                        % 
                        % if ~isempty(alarm_indices)
                        % 
                        %     robustIndex = start_idx + alarm_indices(1) - round((0.025*actualFs - 1)/2);
                        % 
                        %     robustIndex = min(robustIndex, length(selected_Chn));
                        % 
                        %     robustIndex = max(robustIndex, 1); 
                        % 
                        % else
                        % 
                        %     robustIndex = start_idx;
                        % 
                        % end
                        % 
                        % % make the comparison of envelop_hilbert trigger and EMG mean value
                        % % trigger. pick the one comes first.
                        % seg_end = min(start_idx + 5*actualFs - 1, length(seegDiff_smooth));
                        % 
                        % seegSegment = seegDiff_smooth(start_idx:seg_end); % 5s-long smoothed EMG data segment after the trigger signal
                        % 
                        % meanval = mean(seegSegment);
                        % 
                        % t_start = start_idx + 0.25*actualFs;
                        % 
                        % t_end = min(start_idx + 4.5*actualFs - 1, length(seegDiff_smooth));
                        % 
                        % t = t_start:t_end; % 0.25s-4.5s in the segment
                        % 
                        % % 大于1.5倍平均值，说明存在部分trial（trigger>0的t）是有反应的
                        % valid_idx = find(seegDiff_smooth(t) >= 1.5*meanval, 1);
                        % 
                        % if ~isempty(valid_idx)
                        % 
                        %     trigger_pos = min(t(1) - 1 + valid_idx, robustIndex);
                        % 
                        % else
                        % 
                        %     trigger_pos = robustIndex;
                        % 
                        % end
        
                        % 开始构造Gamma3维矩阵
                        
                        
                        A = selected_Chn( ...
                            max(1,index_of_stimulus_onset(trial)-actualFs*extent1): ...
                            min(index_of_stimulus_onset(trial)+actualFs*extent2-1,length(Trig)));
        
                        
        
                        if length(A) == g4
        
                          gamma_epoch(:, trial, chn) = A;
        
                        


                        % if index_of_stimulus_onset(trial) <= Trigger_ind(15)
                        % 
                        %     gamma_epoch_4(:,trial,chn,1) = A;    
                        % 
                        % elseif index_of_stimulus_onset(trial) > Trigger_ind(16) && index_of_stimulus_onset(trial) <= Trigger_ind(30)
                        % 
                        %     gamma_epoch_4(:,trial,chn,2) = A;
                        % 
                        % else 
                        % 
                        %     gamma_epoch_4(:,trial,chn,3) = A;
                        % 
                        % end




                        end
                        
                        
    
    
    
                        % if (floor(((chn-1)*g3+trial)/g3/g2*100)-B==1)
                        % 
                        % count = count+1;
                        % fprintf("%d percent has completed!", count);
                        % 
                        % end
                        % B = floor(((chn-1)*g3+trial)/g3/g2*100);
                    
    
            end
        end
    
        Gamma_epoch_cell{con,cons} = gamma_epoch;
    
        disp(size(gamma_epoch));


    end
end


% 
% for con = 1:length(targetSubjects)
% 
%     con = 
% 
% 
% 
% 
% 
% 
% 
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 


function subInfo = get_subject_info(subjId)
    persistent subjectDB
    if isempty(subjectDB)
        subjectDB = initialize_database();
    end

    idx = find([subjectDB.subjId] == subjId, 1);
    if isempty(idx)
        error('Subject %d not registered', subjId);
    end
    subInfo = subjectDB(idx);
end


function subjectDB = initialize_database()
    GLOBAL_DEFAULTS = struct(...
        'Fs',        2000,...
        'Session_num', [1,2],...
        'Notes',     {{}});
    
    CHANNEL_RULES = struct(...
        'UseChn',  @(x)validateattributes(x,{'numeric'},{'row','positive'}),...
        'EmgChn',  @(x)validateattributes(x,{'numeric'},{'row','numel',2}),...
        'TrigChn', @(x)validateattributes(x,{'numeric'},{'row','numel',3}));

    subjectDB = [
        % Fs=1000
        define_subject(51, 'Fs', 1000, 'TrigChn', [46,47,49], 'UseChn',[1:19,21:37,54:209], 'EmgChn',210:211);
        define_subject(42, 'TrigChn', 32:34, 'UseChn',[1:15,17:31,40:117], 'EmgChn',124:125);
        define_subject(44, 'TrigChn', 32:34, 'UseChn',[1:19,21:35,52:183], 'EmgChn',187:188);
        define_subject(46, 'TrigChn', 46:48, 'UseChn',[1:19,21:37,54:153], 'EmgChn',154:155);
        define_subject(47, 'TrigChn', 46:48, 'UseChn',[1:19,21:37,54:201], 'EmgChn',202:203);
        define_subject(48, 'TrigChn', 39:41, 'UseChn',[1:18,20:30,47:67,80:147,152:161], 'EmgChn',168:169);
        define_subject(49, 'TrigChn', 38:40, 'UseChn',[1:15,17:29,46:181], 'EmgChn',182:183);
        define_subject(50, 'TrigChn', [46,47,49], 'UseChn',[1:19,21:37,54:207], 'EmgChn',208:209);
        
    ];

    [~, idx] = sort([subjectDB.subjId]);
    subjectDB = subjectDB(idx);

    function s = define_subject(id, varargin)
        s = struct('subjId', id);
        for k = 1:2:length(varargin)
            s.(varargin{k}) = varargin{k+1};
        end
        
        s = merge_struct(GLOBAL_DEFAULTS, s);
        
        fields = fieldnames(CHANNEL_RULES);
        for f = 1:length(fields)
            field = fields{f};
            if ~isfield(s, field) || isempty(s.(field))
                error('Subject %d ȱ ٱ Ҫ ֶ : %s', id, field);
            end
            CHANNEL_RULES.(field)(s.(field));
        end
    end
    function S = merge_struct(default, specific)
        S = default;
        fields = fieldnames(specific);
        for i = 1:numel(fields)
            S.(fields{i}) = specific.(fields{i});
        end
    end
end



function alarm = envelop_hilbert_v2(y,Smooth_window,threshold_style,DURATION,gr)
    %% function alarm = envelop_hilbert(y,Smooth_window,threshold_style,DURATION,gr)
    %% Inputs ;
    % y = Raw input signal to be analyzed
    % Smooth_window :this is the window length used for smoothing your signal
    % threshold_style : set it 1 to have an adaptive threshold and set it 0
    % to manually select the threshold from a plot
    % DURATION : Number of the samples that the signal should stay
    % gr = make it 1 if you want a plot and 0 when you dont want a plot
    
    %%%%%%%
    % Tuning parameters for the best results;
    %%%%%%%
    % 1. DURATION is correlated to your sampling frequency, you can use a multiple
    % of your sampling frequency e.g. round(0.050*SamplingFrequency)
    % 2. Smooth_window is correlated to your sampling frequency, you can use a multiple
    % of your sampling frequency e.g. round(0.0500*SamplingFrequency), this is
    % the window length used for smoothing your signal
    
    
    
    
    %% Outputs ;
    % alarm : vector resembeling the active parts of the signal
    %% Method
    % Calculates the analytical signal with the help of hilbert transfrom,
    % takes the envelope and smoothes the signal. Finally , with the help of an
    % adaptive threshold detects the activity of the signal where at least a
    % minimum number of samples with the length of 
    % (DURATION) Samples should stay above the threshold). The threshold is a
    % computation of signal noise and activity level which is updated online.
    
    %% Example and Demo
    % To run demo mode simply execute the following line without any input;
    % Example 1 :
    % alarm = envelop_hilbert()
    % The script generates one artificial signal and analysis that
    % v = repmat([.1*ones(200,1);ones(100,1)],[10 1]); % generate true variance profile
    % y = sqrt(v).*randn(size(v));
    
    % Example 2 : For real world signals with a certain Sampling frequency
    % called (Fs) (In this example a smoothing window with length 200 msec,)
    % alarm = envelop_hilbert(signal,round(0.050*Fs),1,round(0.020*Fs),1)
    
    %% Author : Hooman Sedghamiz 
    % hoose792@student.liu.se
    %(Hooman.sedghamiz@medel.com)
    % Copy right April 2013
    
    % Edited March 2014
    
    %%
    
    % input handling
    if nargin < 5
        gr = 1;
        if nargin < 4
            DURATION = 20; % default
            if nargin < 3
                threshold_style = 1; % default 1 , means it is done automatic
                if nargin < 2
                    Smooth_window = 20; % default for smoothing length
                    if  nargin < 1
                        v = repmat([.1*ones(200,1);ones(100,1)],[10 1]); % generate true variance profile
                        y = sqrt(v).*randn(size(v));
                    end
                end
            end
        end
    end
    
    %% calculate the analytical signal and get the envelope
    test=y(:);
    analytic = hilbert(test);
    env = abs(analytic);
    
    %% take the moving average of analytical signal
    %env=movingav(env,70,0);
    env = conv(env, ones(1,Smooth_window)/Smooth_window);%smooth
    env = env(:) - mean(env); % get rid of offset
    env = env/max(env); %normalize
    
    %% threshold the signal
    if threshold_style == 0
        hg=figure;plot(env);title('Select a threshold on the graph')
        [~, THR_SIG] =ginput(1);
        close(hg);
    end
    %DURATION = 20;
    
    h=1;
    alarm =zeros(1,length(env));
    if threshold_style
        THR_SIG = 4*mean(env);
    end
    nois = mean(env)*(1/3); % noise level
    threshold = mean(env); % signal level
    
    thres_buf  = [];
    nois_buf = [];
    THR_buf = zeros(1,length(env));
    
    for i = 1:length(env)-DURATION
      if env(i:i+DURATION) > THR_SIG 
          % alarmx(h) = i;
          % alarmy(h) = env(i);
          alarm(i) = max(env);
          threshold = 0.2*mean(env(i:i+DURATION)); % update threshold 10% of the maximum peaks found
          h = h + 1;
      else
          if mean(env(i:i+DURATION)) < THR_SIG
          nois = mean(env(i:i+DURATION)); %update noise
          else
              if ~isempty(nois_buf)
                  nois = mean(nois_buf);
              end
          end
      end 
      
      thres_buf = [thres_buf, threshold];
      nois_buf = [nois_buf, nois];
      
      if h > 1
      THR_SIG = nois + 0.50*(abs(threshold - nois)); %update threshold
      end
      THR_buf(i) = THR_SIG;
    end 
    
    if gr
    figure,ax(1)=subplot(211);plot(test/max(test)),hold on,plot(alarm/(max(alarm)),'r','LineWidth',2.5),
    hold on,plot(THR_buf,'--g','LineWidth',2.5);
    title('Raw Signal and detected Onsets of activity');
    legend('Raw Signal','Detected Activity in Signal','Adaptive Treshold',...
        'orientation','horizontal');
    grid on;axis tight;
    ax(2)=subplot(212);plot(env);
    hold on,plot(THR_buf,'--g','LineWidth',2.5),
    hold on,plot(thres_buf,'--r','LineWidth',2),
    hold on,plot(nois_buf,'--k','LineWidth',2),
    title('Smoothed Envelope of the signal(Hilbert Transform)');
    legend('Smoothed Envelope of the signal(Hilbert Transform)','Adaptive Treshold',...
        'Activity level','Noise Level','orientation','horizontal');
    linkaxes(ax,'x');
    zoom on;
    axis tight;
    grid on;
    end

end



