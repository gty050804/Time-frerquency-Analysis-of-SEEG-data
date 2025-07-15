% 
con=0;

g1 = size(targetSubjects,2);

reactive_locations = cell(g1,2);

non_reactive_locations = cell(g1,2);

r_obs_cell = cell(g1,2);

for subjId = targetSubjects

    fprintf("Identifying location of P%d \n",subjId);

    subInfo = get_subject_info(subjId);

    Fs = subInfo.Fs;   actualFs = Fs;

    con = con+1;

    sessionNo = subInfo.Session_num;

    for cons = sessionNo

    gamma_interest = Gamma_epoch_cell{con,cons};



    % extent1表示的是Trigger前的时间，extent2表示的是Trigger后的时间。
 
    baseline_low_limit = 0.8;

    baseline_high_limit = 1;          % 设置基线时间选取范围（这里即为3s*0.25~0.75=0.75~2.25s的1500ms的时间范围内）

    baseline = gamma_interest(round(extent1*baseline_low_limit*actualFs) ...
        :round(extent1*baseline_high_limit*actualFs),:,:);

    task_low_limit = 0;

    task_high_limit = 0.1;           % 设置任务相关时间选取范围（这里即为6s*0~0.25=0~1.5s的1500ms的时间范围）

    task = gamma_interest(round(extent1*actualFs+extent2*task_low_limit*actualFs)+1 ...
        :round(extent1*actualFs+extent2*task_high_limit*actualFs),:,:);

    % 基线与任务相关的样本数量

    N_samples_of_baseline = size(baseline,1);

    N_samples_of_task = size(task,1);


    %% permutation テスト  （排列测试）

    number_of_repetitions = 1000;   % 实验重复次数

    reverseStr = '';

    number_of_locations = size(gamma_interest,3);

    r_obs = zeros(1,number_of_locations);

    r_pdf = zeros(number_of_repetitions,number_of_locations);

    p_value = zeros(1,number_of_locations);

    chiver = zeros(1,number_of_locations);

    for chn = 1:number_of_locations

        % displaying message
        msg = sprintf('\n Processing channel %d/%d', chn, number_of_locations);
        fprintf([reverseStr, msg]);
        

        % 检验各有效Trigger的基线和任务相关阶段均值相关性
        % calculating distributions x (mean of each baseline epoch) and y (mean of each task epoch)
        x = mean(baseline(:, :, chn), 1); y = mean(task(:, :, chn), 1);
        N = length(x);
        M = length(y);
        labels = [ones(1,N)*-1, ones(1,M)]; % creates a vector of labels (-1 for baseline, +1 for task)
        r_obs(chn) = corr([x y]', labels', 'type','spearman'); % computes r_obs

        % calculating the distribution
        % r_pdf（r_pdf：相关性不大的，标签啥的都是瞎配的；r_obs：按理来说有相关性的，标签都是对应的）
        tmp = [x y];
        for i = 1:number_of_repetitions
            
            w = tmp(randperm(size(tmp, 2))); % permutating the means of baseline and task epochs
            r_pdf(i, chn) = corr(w', labels', 'type','spearman'); % computing correlation coefficient (i.e., one sample of the distribution r_pdf)
            
        end

        p_value(chn) = 2*normcdf(-abs(r_obs(chn)), mean(r_pdf(:, chn), 1), std(r_pdf(:, chn), 0, 1));
        chiver(chn) = chi2gof(r_pdf(:,chn)); % reture value=0 means the data match the normality, otherwise, doesn't match the normality.

    end

    % 已知信息：r_obs： 目标相关性
    %          r_pdf： 随机测试相关性
    %          p_value：p值检验（p越小，相关程度越强）
    %          chiver：随机测试正态分布检验

    %% selecting locations with p-values smaller than 0.05 after Bonferroni-correcting for the number of tests (i.e., number of locations)
    % To do Bonferroni correction, we divide 0.05 (significance level) by
    % the number of locations.
    % reactive_locations = find(p_value < 0.05); 
    % reactive_locations = find(p_value < 0.05/number_of_locations);
    [~,reactive_locations_1] = min(p_value);

    %% plots r_pdf and r_obs for one reactive location and one non-reactive location

    figure
    histogram(r_pdf(:, reactive_locations_1)),
    hold on
    x1 = abs(r_obs(reactive_locations_1));
    y1 = get(gca,'ylim');
    plot([x1 x1],y1, 'r')
    title('Distribution of correlation coefficients observed by chance (r pdf, in blue) versus actual observed coefficient (r obs, red line) for a reactive location')
    
    
    % non_reactive_locations = setdiff(1:number_of_locations, reactive_locations); % determining non-reactive locations
    [~,non_reactive_locations_1] = max(p_value);


    figure
    histogram(r_pdf(:, non_reactive_locations_1)),
    hold on
    x1 = abs(r_obs(non_reactive_locations_1));
    y1 = get(gca,'ylim');
    plot([x1 x1],y1, 'r')
    title('Distribution of correlation coefficients observed by chance (r pdf, in blue) versus actual observed coefficient (r obs, red line) for a non-reactive location')

    reactive_locations{con,cons} = reactive_locations_1;
    non_reactive_locations{con,cons} = non_reactive_locations_1;
    r_obs_cell{con,cons} = r_obs;
    end


end



























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
        'TrigChn', @(x)validateattributes(x,{'numeric'},{'row','numel',5}));

    subjectDB = [
        % Fs=1000
        define_subject(2,  'Fs',1000, 'UseChn',[1:19,21:37,43:44,47:129], 'EmgChn',145:146, 'TrigChn',38:42)
        define_subject(3,  'Fs',1000, 'Session_num',[1,3], 'UseChn',[1:19,21:37,44:45,48:189], 'EmgChn',192:193, 'TrigChn',38:42)
        define_subject(4,  'Fs',1000, 'Session_num',[2,3], 'UseChn',[1:19,21:37,43:44,47:68], 'EmgChn',75:76, 'TrigChn',38:42)
        define_subject(5,  'Fs',1000, 'Session_num',[1,3], 'UseChn',[1:19,21:37,43:44,47:150,151:166,167:186],...
                         'EmgChn',187:188, 'TrigChn',38:42, 'Notes',{{'N1 has some drifts during the entire session, and need to be removed'}})
        define_subject(7,  'Fs',1000, 'UseChn',[1:19,21:37,44:45,48:126,128:153], 'EmgChn',162:163, 'TrigChn',38:42)
        define_subject(8,  'Fs',1000, 'UseChn',[1:36,53:141,143:186], 'EmgChn',193:194, 'TrigChn',38:42,...
                         'Notes',{{'N1 has some drifts during the entire session, and need to be removed'}})
        define_subject(9,  'Fs',1000, 'UseChn',[1:19,21:37,44:45,48:123], 'EmgChn',124:125, 'TrigChn',38:42)
        define_subject(20, 'Fs',1000, 'UseChn',[1:19,21:37,47:48,51:120], 'EmgChn',127:128, 'TrigChn',39:43)
        define_subject(21, 'Fs',1000, 'UseChn',[1:19,21:37,46:47,50:129], 'EmgChn',136:137, 'TrigChn',38:42)
        define_subject(35, 'Fs',1000, 'UseChn',[1:19,21:31,40:41,44:57,60:87,90:151],...
                         'EmgChn',88:89, 'TrigChn',32:36)

        % Fs=2000
        define_subject(10, 'TrigChn',34:38, 'UseChn',[1:19,21:33,39:60,63:216], 'EmgChn',61:62)
        define_subject(13, 'UseChn',[1:19,21:35,52:119], 'EmgChn',126:127, 'TrigChn',44:48)
        define_subject(14, 'UseChn',[1:17,19:35,44:139], 'EmgChn',142:143, 'TrigChn',36:40,...
                         'Notes',{{'not the correct EMG'}})
        define_subject(16, 'UseChn',[1:19,21:37,46:179], 'EmgChn',186:187, 'TrigChn',38:42)
        define_subject(17, 'UseChn',[1:19,21:37,46:153], 'EmgChn',160:161, 'TrigChn',38:42)
        define_subject(18, 'UseChn',[1:19,21:35,52:161], 'EmgChn',182:183, 'TrigChn',44:48)
        define_subject(19, 'UseChn',[1:19,21:37,46:146], 'EmgChn',153:154, 'TrigChn',38:42)
        define_subject(22, 'UseChn',[1:17,19:33,42:159], 'EmgChn',160:161, 'TrigChn',34:38)
        define_subject(23, 'UseChn',[1:16,18:34,43:161,168:213], 'EmgChn',214:215, 'TrigChn',35:39)
        define_subject(24, 'UseChn',[1:14,16:30,47:145], 'EmgChn',146:147, 'TrigChn',39:43,...
                         'Notes',{{'miss K5 & A3 ,add two virtual channels for P24, and should be defined as bad channels'}})
        define_subject(25, 'UseChn',[1:15,17:30,40:146], 'EmgChn',147:148, 'TrigChn',32:36,...
                         'Notes',{{'H9 is missing, create one virtual channel for P25'}})
        define_subject(26, 'UseChn',[1:17,19:33,43:163], 'EmgChn',164:165, 'TrigChn',35:39,...
                         'Notes',{{'L6 is missing for P26'}})
        define_subject(29, 'UseChn',[1:15,17:29,38:119], 'EmgChn',120:121, 'TrigChn',30:34)
        define_subject(30, 'UseChn',[1:15,17:29,38:103,106:119], 'EmgChn',120:121, 'TrigChn',30:34)
        define_subject(31, 'UseChn',[1:17,19:33,42:81], 'EmgChn',82:83, 'TrigChn',34:38)
        define_subject(32, 'UseChn',[1:19,21:37,46:47,50:67], 'EmgChn',68:69, 'TrigChn',38:42)
        define_subject(34, 'UseChn',[1:15,17:31,43:114], 'EmgChn',117:118, 'TrigChn',35:39)
        define_subject(36, 'UseChn',[1:15,17:25,34:126], 'EmgChn',127:128, 'TrigChn',26:30)
        define_subject(37, 'UseChn',[1:15,17:23,32:73], 'EmgChn',74:75, 'TrigChn',24:28)
        % define_subject(38, 'UseChn',[1:4,7:19,21:37,46:253],
        % 'EmgChn',254:255, 'TrigChn',38:42) ECoG
        define_subject(39, 'UseChn',[1:19,21:35,44:135], 'EmgChn',138:139, 'TrigChn',36:40)
        define_subject(41, 'UseChn',[1:19,21:37,54:207], 'EmgChn',210:211, 'TrigChn',46:50)
        % define_subject(43, 'UseChn',[1:19,21:37,54:225],
        % 'EmgChn',228:229, 'TrigChn',46:50) ECoG
        define_subject(45, 'UseChn',[1:19,21:37,46:181], 'EmgChn',182:183, 'TrigChn',38:42)

        % Fs=500
        define_subject(11, 'Fs',500, 'UseChn',[1:19,21:35,44:45,48:143,146:151,154:209],...
                         'EmgChn',216:217, 'TrigChn',36:40)
        define_subject(12, 'Fs',500, 'UseChn',[1:19,21:37,43:44,47:102], 'EmgChn',109:110, 'TrigChn',38:42)
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











