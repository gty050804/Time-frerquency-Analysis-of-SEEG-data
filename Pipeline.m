% Copyright (c) 2025 Shanghai Jiao Tong University. All rights reserved.
% Created: 2025.3.21
% Last Modified: 2025.7.5

% A_M:
% 目前5通道可用P：7 8 9 10 12 13 16 17 18 19(需要改名字）20 21 22 29 30 31 36 38(ECog) 45 
% 目前5通道不可用P：14 15 34(Only 10/90 left) 
% 3通道P：42 44 46 47 48 49 50 51
% 未知P：6 23 

% V_M:
% 目前5通道可用P：7 8 9 10 11 12 16 17 18 19 20 21? 22 29? 34?? 43 45? 
% 目前5通道不可用P：14 15 
% 3通道P：
% 未知P：6 23 36(S1出现100个trigger） 
% 1Session: 13 

clc,clear

tic
targetSubjects = 19;

AV = 1;   % 实验选项

% P7 肌电图噪声巨大
% P9 部分通道效果较好
% P10 




% P2 3 4 5数据在文件中没有出现，并且SessionNum出现一些不同，此处忽略。
% 因此默认所有P的SessionNum都是1:2。

% emg分析时，若emg_trig（最后一列）和trigger（倒数第二列）的相同行差值小于3*actualFs，
% 则说明在想想运动的过程中手臂出现了运动，对于这部分Trigger需要予以剔除。


count=0;

[~, name] = system('hostname');
% if strcmp(strip(name), 'ACE1999')
    config.root_dir = fullfile('D:', 'Project','Code_Tutorial','03_Imaginary_Gesture_Coding','Coding_Gestures_Imaginary');
    config.raw_dir = fullfile(config.root_dir, '1_Raw_Data_All');
    config.result_dir = fullfile(config.root_dir, 'preprocessed_data');
    config.EleCTX_dir = fullfile(config.root_dir, '3_Brain_Electrodes');
    config.EMG_dir = fullfile(config.root_dir,"EMG_data");
    safe_mkdir(config.result_dir)
% else
    % error('Please make a config for yourself')
% Example of your dir setting
% root_dir
% │  
% ├─EleCTX_dir
% │  ├─P10
% │  │      WholeCortex.mat
% │  │      SignalChanel_Electrode_Registration.mat
% │  │      electrode_raw.mat
% │  │      electrodes_Final_Norm.mat
% │  │      电极定位图.jpg
% │  │      
% │  ├─P11
% │  │      WholeCortex.mat
% │  │      SignalChanel_Electrode_Registration.mat
% │  │      electrode_raw.mat
% │  │      electrodes_Final_Norm.mat
% │  │      电极定位图.jpg
% │   
% │  
% ├─raw_dir       
% │  ├─P39
% │  │  │  Configuration_Setup.m
% │  │  │  
% │  │  └─1_Raw_Data_Transfer
% │  │          Save_Raw_SEEG_Data.m
% │  │          P39_H1_1_Raw.mat
% │  │          P39_H1_2_Raw.mat
% │  │          
% │  ├─P34
% │  │  │  Configuration_Setup.m
% │  │  │  
% │  │  └─1_Raw_Data_Transfer
% │  │          Save_Raw_SEEG_Data.m
% │  │          P34_H1_1_Raw.mat
% │  │          P34_H1_2_Raw.mat
% │ 
% │              
% ├─result_dir
% │  ├─P5
% │  │      preprocessed1.mat
% │  │      preprocessed2.mat
% │  │      
% │  ├─P7
% │  │      preprocessed1.mat
% │  │      preprocessed2.mat
% end

%%%%%%%%%%%%%%%%% Select your target subject id for preprocessing  %%%%%%%%%%%%%%%
% targetSubjects = [2:4, 10, 13, 17, 25, 29, 32, 41, 45];
% targetSubjects = [5, 7:9, 11:12, 14, 16, 18:24, 26, 30:31, 34:37, 39];
sessionNum = 2;

Subjectcell = cell(length(targetSubjects),sessionNum);

Trigger_ind_cell = cell(length(targetSubjects),sessionNum);

% process pipline.
for subjId = targetSubjects
    fprintf('\nProcessing subject %d...', subjId);

    count = count+1;



    if ~ismember(subjId,[42,44,46,47,48,49,50,51])
        subInfo = get_subject_info(subjId);
    else
        subInfo = get_subject_info_3(subjId);
    end

    sessionNum = numel(subInfo.Session_num);

    cycle = 0;
    
    while true

    % downsample SEEG & EMG, select good SEEG channels, process EMG, and detect triggers
        [Datacell, good_channels, actualFs] = preprocess_stage1(config, subjId, subInfo,AV);

        cycle = cycle+1;
    
        sum_1 = Datacell{1,1}(:,end);
    
        sum_2 = Datacell{1,2}(:,end);
    
        if sum(sum_1) == 45 && sum(sum_2) == 45
            break;
        end

        fprintf("Failure in complete trigger detection. Processing trial:%d",cycle+1);


    end



    % disp(size(Datacell));

    % 由于P10出现上升沿检测鲁棒性低的问题，此处采用手动上升沿检测。
    
    % if subjId == 10
    % 
    %     Trig_1 = zeros(size(Datacell{count,1},1),1);
    % 
    %     % disp(size(Trig_1));
    % 
    %     Trig_2 = zeros(size(Datacell{count,2},1),1);
    % 
    %     Trig_idx_1 = [27299,36463,45628,54826,64057,73322,82569,91801,101065,110313,119611,128909,138240,147571,156943,166200,175598,185196,194527,203975,213372,222737,232168,241616,251080,260495,269926,279407,288838,298303,307834,317332,326846,336427,345908,355423,365021,374518,384649,394330,403978,413543,423140,432738,442353];
    % 
    %     Trig_idx_2 = [21896,31034,40215,49446,58694,67942,77173,86404,95669,104933,114231,123512,132793,142058,151372,160720,170068,179416,188797,198162,207509,216891,226272,235653,245184,254598,264196,273644,283059,292523,302004,311569,321050,330614,340112,349610,359174,368689,378503,388034,397632,407230,416827,426409,436190];
    % 
    %     Trig_1(Trig_idx_1,1) = 1;
    % 
    %     Trig_2(Trig_idx_2,1) = 1;
    % 
    %     Datacell{count,1}(:,end) = Trig_1';
    % 
    %     Datacell{count,2}(:,end) = Trig_2';
    % 
    % 
    % end


    % process SEEG (bandpass, notch, re-reference), and modify triggers using EMG
    [processed_cell,Trigger_ind] = preprocess_stage2(config, subjId, Datacell, good_channels, actualFs);
    
    
    for n = 1:sessionNum
        Subjectcell{count,n} = processed_cell{n};
        Trigger_ind_cell{count,n} = Trigger_ind{n};
    end

end

toc

EleCTX_dir = fullfile(config.EleCTX_dir, sprintf('P%d', subjId));
load(fullfile(EleCTX_dir, 'electrodes_Final_Norm.mat'), 'elec_Info_Final_wm');
name = elec_Info_Final_wm.name;

goodname = name(good_channels);













%% 
function [Datacell, good_channels, actualFs] = preprocess_stage1(config, subjId, subInfo,AV)
    fprintf('\n-- Stage 1 --');
    
    Fs = subInfo.Fs;
    sessionNum = numel(subInfo.Session_num);
    Datacell = cell(1, sessionNum);
    goodChsCache = cell(1, sessionNum);

    for sessionIdx = 1:sessionNum
        % load session data
        fprintf(' Session %d/%d', sessionIdx, sessionNum);

        if AV == 2

            if subjId == 19 
    
                dataFile = fullfile(config.raw_dir, sprintf('P%d', subjId), '1_Raw_Data_Transfer', ...
                sprintf('P%d_V_M_%d.mat', subjId, subInfo.Session_num(sessionIdx)));
    
            else
    
                dataFile = fullfile(config.raw_dir, sprintf('P%d', subjId), '1_Raw_Data_Transfer', ...
                sprintf('P%d_V_M_%d_Raw.mat', subjId, subInfo.Session_num(sessionIdx)));
    
            end

        elseif AV == 1

            if subjId == 19 
    
                dataFile = fullfile(config.raw_dir, sprintf('P%d', subjId), '1_Raw_Data_Transfer', ...
                sprintf('P%d_A_M_%d.mat', subjId, subInfo.Session_num(sessionIdx)));
    
            else
    
                dataFile = fullfile(config.raw_dir, sprintf('P%d', subjId), '1_Raw_Data_Transfer', ...
                sprintf('P%d_A_M_%d_Raw.mat', subjId, subInfo.Session_num(sessionIdx)));
    
            end


        end
        if ~exist(dataFile, 'file')
            error('Data file not found: %s', dataFile);
        end
        load(dataFile, 'Data');

        % resampling
        [data, actualFs] = handle_resampling(Data, Fs, [subInfo.UseChn, subInfo.EmgChn]);
    
        % select good channals
        goodChsCache{sessionIdx} = select_good_channels(data(:, subInfo.UseChn), actualFs, 10);

        % process EMG
        [emgProc, emgDiff] = process_emg(data(:, subInfo.EmgChn), actualFs);
        
        % detect triggers 
        trigger_labels = detect_triggers(data(:, subInfo.TrigChn),actualFs, subjId,sessionIdx,AV);

        % disp(size(data(:, subInfo.UseChn)));
        % disp(size(emgProc));
        % disp(size(emgDiff));
        % disp(size(trigger_labels));
        
        % (seeg x N, emg x 2, emgdiff x 1, feature x 1)
        Datacell{sessionIdx} = [data(:, subInfo.UseChn), emgProc, emgDiff, trigger_labels];   
    end
    
    good_channels = union(goodChsCache{:});
end

function [output,Trigger_cell_cell] = preprocess_stage2(config, subjId, Datacell, good_channels, actualFs)
    fprintf('\n-- Stage 2 --\n');
    Trigger_cell_cell = cell(2,1);
    sessionNum = numel(Datacell);
    emg_trig_cell = cell(1,2);
    ranking = zeros(1,90);
    emg_cell = cell(2,45);
    for sessionIdx = 1:sessionNum
        fprintf(' Session %d/%d:', sessionIdx, sessionNum);
        sessionData = Datacell{sessionIdx};  % seeg*n + emg*2 + emgdiff + trigger_exe
        trigger_idx = sessionData(:, end);
        EMG = sessionData(:, end-3:end-2);
        SEEG = sessionData(:, 1:end-4);

        
        
        % bandpass and notch
        k = 2;
        OME = [0.5, 400];
        SEEG = cFilterD_EEG(SEEG, actualFs, k, OME); % with a Notch filter
        
        if subjId ~= 38 && subjId ~=43
        % re-reference
            [SEEG_referenced, orderedChannelIdx] = Laplacian_reRef(config, SEEG, good_channels, subjId);
        else

            [SEEG_referenced, orderedChannelIdx] = CAR_reRef(config, SEEG, good_channels, subjId);

        end
        % SEEG channel re-order
        SEEG_referenced = SEEG_referenced(:, orderedChannelIdx);
        [~, good_channels] = ismember(good_channels, orderedChannelIdx);

        % alignment of feature_lable and EMG.
        emgDiff = sessionData(:, end-1); % emg diff data


        % 滤波器处理emgDiff
        fc = 50; % 高通截止频率（Hz）
        fs = actualFs; % 采样频率（Hz）

        % Butterworth高通滤波器
        [b, a] = butter(4, fc/(fs/2), 'high');

        % 应用滤波器（零相位滤波）
        emgDiff = filtfilt(b, a, emgDiff);


        % emgDiff = process_emg_2(emgDiff,actualFs);



        % plot((1:length(emgDiff))/actualFs, emgDiff');



        emgDiff_smooth = smooth(abs(emgDiff), 0.025*actualFs);


        % plot((1:length(emgDiff_smooth))/actualFs, emgDiff_smooth');
        EMG_trigger = zeros(size(sessionData,1), 1);
        trigger = find(trigger_idx); % search for the trigger position and label
        Trigger_cell_cell{sessionIdx,1} = trigger;





        
        for trial = 1:length(trigger)  % the segment number (i th)
            start_idx = trigger(trial);
            end_idx = min(start_idx + 5*actualFs, length(emgDiff));
            emgSegment = emgDiff(start_idx:end_idx); % 5s-long EMG data segment after the trigger signal

            % emgSegment = gpuArray(emgSegment);

            % detect event
            alarm = envelop_hilbert_v2(emgSegment, round(0.025*actualFs), 1, round(0.05*actualFs), 0);

            % alarm = gather(alarm);
            alarm_indices = find(alarm == 1);
            
            if ~isempty(alarm_indices)

                robustIndex = start_idx + alarm_indices(1) - round((0.025*actualFs - 1)/2);

                robustIndex = min(robustIndex, length(emgDiff));

                robustIndex = max(robustIndex, 1); 
                
            else
                robustIndex = start_idx;
                
            end
            
            % make the comparison of envelop_hilbert trigger and EMG mean value
            % trigger. pick the one comes first.
            seg_end = min(start_idx + 5*actualFs - 1, length(emgDiff_smooth));
            emgSegment = emgDiff_smooth(start_idx:seg_end); % 5s-long smoothed EMG data segment after the trigger signal
            meanval = mean(emgSegment);

            t_start = start_idx + 0.25*actualFs;
            t_end = min(start_idx + 4.5*actualFs - 1, length(emgDiff_smooth));
            t = t_start:t_end; % 0.25s-4.5s in the segment
            valid_idx = find(emgDiff_smooth(t) >= 10*abs(meanval), 1);
            if ~isempty(valid_idx)
                trigger_pos = min(t(1) - 1 + valid_idx, robustIndex);
            else
                trigger_pos = robustIndex;
      
            end
            EMG_trigger(trigger_pos) = trigger_idx(trigger(trial));
        end

        % (eegdata, 2*EMG, 1*fealabel, 1*EMG_trigger)
        Datacell{sessionIdx} = [SEEG_referenced, EMG, trigger_idx, EMG_trigger]; 

        temp = Datacell{sessionIdx}(:,end-1);

        EMG_trig = find(EMG_trigger>0);

        

        % disp(length(EMG_trig));



        t_diff = EMG_trig - trigger;

        emg_trig_cell{1,sessionIdx} = t_diff;

        disp(size(t_diff));

        to_be_deleted = trigger(t_diff<3*actualFs);

        % disp(length(to_be_deleted));

        disp(length(trigger));

        fprintf("%d triggers have been eliminated due to EMG detection \n",length(to_be_deleted));

        temp(to_be_deleted) = 0;


        % emg = find(EMG_trigger>0);
        % to_be_deleted = find(emg - trigger<3*actualFs);
        % to_be_deleted = trigger_idx(to_be_deleted);
        % temp(to_be_deleted) = 0;
        Datacell{sessionIdx}(:,end-1) = temp;

        Datacell{sessionIdx} = Datacell{sessionIdx}(:,1:end-1);

        Datacell{sessionIdx} = Datacell{sessionIdx}(:,[good_channels,end-2,end-1,end]);


        

        % ranking = zeros(1,90);




        

        % for k = 1:3
        % 
        %     index = [1:15,46:60] + (k-1)*15;
        % 
        %     to_be_rank = EMG_trig(index);
        % 
        %     [~,~,rank] = unique(to_be_rank);
        % 
        %     for j = 1:length(rank)
        % 
        %         ranking((ceil(j/15)-1)*45+mod(j,15)) = rank(j);
        % 
        %     end
        % 
        % 
        % 
        % end

        set(0,'DefaultFigureVisible', 'on');
        
        if sessionIdx == 1

            figure

        end

        for ind = 1:length(trigger)

            if ind ~= length(trigger)

                % disp(min(trigger(ind)-3*actualFs:trigger(ind)+6*actualFs));
                % 
                % disp(max(trigger(ind)-3*actualFs:trigger(ind)+6*actualFs));
                if trigger(ind)-3*actualFs <= 0

                    continue;

                end

                emg_cell{sessionIdx,ind} = emgDiff(trigger(ind)-3*actualFs:trigger(ind)+6*actualFs);

            else

                emg_cell{sessionIdx,ind} = emgDiff(trigger(ind)-3*actualFs:min(trigger(ind)+6*actualFs,length(emgDiff)));

            end

            emg_cell{sessionIdx,ind} = (emg_cell{sessionIdx,ind}-mean(emg_cell{sessionIdx,ind}))./sqrt(var(emg_cell{sessionIdx,ind},1));

            emg_cell{sessionIdx,ind} = emg_cell{sessionIdx,ind}/20;

            

        end



    end



    % 现在有的数据：1、emg_cell   2行45列   单cell存储1trigger数据   长度9000的行向量
                  %2、emg_trig_cell   1行2列  待排序数据  单cell存储45个trigger的emg响应时间
     for ka = 1:3

            % disp("Complete");

            id = (ka-1)*15+(1:15);

            to_be_rank = [emg_trig_cell{1,1}(id);emg_trig_cell{1,2}(id)];

            % disp(length(emg_trig_cell{1,1}));

            [~,~,rank] = unique(to_be_rank);

            for m = 1:30

                idx = find(rank == m);

                if length(idx) > 1

                    rank(rank>m) = rank(rank>m)+length(idx)-1;

                    for n = 2:length(idx)

                        rank(idx(n)) = rank(idx(n)) + n-1;

                    end

                end

            end

            % for j = 1:length(rank)
            % 
            %     disp((ceil(j/15)-1)*45+mod(j,15));
            % 
            %     ranking((ceil(j/15)-1)*45+mod(j,15)) = rank(j);
            % 
            % end

            for j = 1:length(rank)

                % disp(j+(ka-1)*15+(ceil(j/15)-1)*30);
                % 
                % disp(rank(j));

                ranking(j+(ka-1)*15+(ceil(j/15)-1)*30) = rank(j);

            end

           

            


     end

     % disp(ranking);

     for j = 1:90

         to_be_transformed = ranking(j) + (mod(ceil(j/15)-1,3))*30 + ...
                                     emg_cell{ceil(j/45),mod(j-1,45)+1};

         len = round(length(to_be_transformed)*2/9);

         shifted = to_be_transformed(1:len);

         main = to_be_transformed(len+1:end);

         merge = [main;shifted];

         % disp(emg_trig_cell{1,ceil(j/45)}(mod(j-1,45)+1))
         % 
         % disp(ranking(j)+(mod(ceil(j/15)-1,3))*30)

         % set(0,'DefaultFigureVisible', 'on');

         plot((1:length(emg_cell{1,ind}))/actualFs,merge,'k');

         hold on

         scatter((emg_trig_cell{1,ceil(j/45)}(mod(j-1,45)+1))/actualFs+1,ranking(j)+(mod(ceil(j/15)-1,3))*30,5,'red','filled');
         % 
         % hold on

         if j == 90

                axis([-0.5,9.5,-5,95]);

                hold on

         end




     end

     x_loc1 = round(length(to_be_transformed)/9)/actualFs;

     x_loc2 = round(length(to_be_transformed)/9*4)/actualFs;

     x_loc3 = round(length(to_be_transformed)/9*7)/actualFs;

     y_min = -5;

     y_max = 95;

     line([x_loc1, x_loc1], [y_min, y_max], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1);

     hold on

     line([x_loc2, x_loc2], [y_min, y_max], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1);

     hold on

     line([x_loc3, x_loc3], [y_min, y_max], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1);

     hold on

     title("EMG records");

     xlabel("time(s)");

     ylabel("trial");

     save_EMG = config.EMG_dir;

     filename = sprintf("P%d",subjId);
    
     if ~exist(save_EMG, 'dir')
         mkdir(save_EMG);
     end

     saveas(gcf, fullfile(save_EMG, filename),'png');


    output = Datacell;




    % save data file.
    save_path = fullfile(config.result_dir, sprintf('preprocessed_P%d.mat', subjId));
    save(save_path, 'Datacell', 'good_channels', 'actualFs', '-v7.3');
end

%%
function safe_mkdir(dirPath)
    if ~exist(dirPath, 'dir')
        mkdir(dirPath);
    end
end

function subjectDB_3 = initialize_database_3()
    GLOBAL_DEFAULTS = struct(...
        'Fs',        2000,...
        'Session_num', [1,2],...
        'Notes',     {{}});
    
    CHANNEL_RULES = struct(...
        'UseChn',  @(x)validateattributes(x,{'numeric'},{'row','positive'}),...
        'EmgChn',  @(x)validateattributes(x,{'numeric'},{'row','numel',2}),...
        'TrigChn', @(x)validateattributes(x,{'numeric'},{'row','numel',3}));

    subjectDB_3 = [
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

    [~, idx] = sort([subjectDB_3.subjId]);
    subjectDB_3 = subjectDB_3(idx);

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
        % define_subject(51,'Fs',1000, 'UseChn',[1:19,21:37,54:209], 'EmgChn',210:211, 'TrigChn',[46,47,49])

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
        define_subject(38, 'UseChn',[1:4,7:19,21:37,46:253],'EmgChn',254:255, 'TrigChn',38:42) 
        define_subject(39, 'UseChn',[1:19,21:35,44:135], 'EmgChn',138:139, 'TrigChn',36:40)
        define_subject(41, 'UseChn',[1:19,21:37,54:207], 'EmgChn',210:211, 'TrigChn',46:50)
        define_subject(43, 'UseChn',[1:19,21:37,54:225],'EmgChn',228:229, 'TrigChn',46:50) 
        define_subject(44, 'UseChn',[1:19,21:35,52:183],'EmgChn',187:188, 'TrigChn',44:48) 
        define_subject(45, 'UseChn',[1:19,21:37,46:181], 'EmgChn',182:183, 'TrigChn',38:42)
        % define_subject(46, 'UseChn',[1:19,21:37,54:153], 'EmgChn',154:155, 'TrigChn',46:48)
        % define_subject(47, 'UseChn',[1:19,21:37,54:201], 'EmgChn',202:203, 'TrigChn',46:48)
        % define_subject(48, 'UseChn',[1:18,20:30,47:67,80:147,152:161], 'EmgChn',168:169, 'TrigChn',39:41)
        % define_subject(49, 'UseChn',[1:15,17:29,46:181], 'EmgChn',182:183, 'TrigChn',38:40)
        % degine_subject(50, 'UseChn',[1:19,21:37,54:207], 'EmgChn',208:209, 'TrigChn',[46,47,49])
        % 

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

function subInfo = get_subject_info_3(subjId)
    persistent subjectDB_3
    if isempty(subjectDB_3)
        subjectDB_3 = initialize_database_3();
    end

    idx = find([subjectDB_3.subjId] == subjId, 1);
    if isempty(idx)
        error('Subject %d not registered', subjId);
    end
    subInfo = subjectDB_3(idx);
end

function [data, actualFs] = handle_resampling(rawData, baseFs, channel)
    persistent b a baseFs_cache

    data = double(rawData');
    if baseFs == 2000
        if isempty(baseFs_cache)
            cutoff = 500;
            [b, a] = butter(4, cutoff/(2000/2), 'low');
            baseFs_cache = 2000;
        end


        data(:, channel) = filtfilt(b, a, data(:, channel));
        data = data(1:2:end, :);
        actualFs = 1000;
    else
        actualFs = baseFs;
    end
end

function good_chs = select_good_channels(signal, sampling_rate, para)
    % gets the 50 Hz noise power using an IIR peak filter
    [b, a] = iirpeak(50/(sampling_rate/2), 0.001);
    noise_level = mean(sqrt(filter(b,a, signal).^2), 1);

    % the good channels are defined as those which power in the line noise is
    % not significantly different from the median line noise across all 
    % channels. We arbitrarily define significance as the median of the 
    % noise across channels + 10 times its median absolute deviation: 
    idx = noise_level < (median(noise_level) + para*mad(noise_level, 1));
    good_chs = find(idx);
end

function [emgProc, emgDiff] = process_emg(emgRaw, Fs)
    persistent bComb aComb bBand aBand Fs_cache
    
    if isempty(Fs_cache) || Fs ~= Fs_cache
        F0 = 50; q = 30;
        n = round(Fs/F0);
        bw = (F0/(Fs/2))/q;
        [bComb, aComb] = iircomb(n, bw, 'notch');
        bpRange = [1.5, 150];
        [bBand, aBand] = butter(4, bpRange/(Fs/2));
        Fs_cache = Fs;
    end  
    
    emgProc = filtfilt(bComb, aComb, emgRaw);
    emgProc = filtfilt(bBand, aBand, emgProc);
    emgDiff = emgProc(:,1) - emgProc(:,2);
end

% function emgProc = process_emg_2(emgRaw, Fs)
%     persistent bComb aComb bBand aBand Fs_cache
% 
%     if isempty(Fs_cache) || Fs ~= Fs_cache
%         F0 = 50; q = 30;
%         n = round(Fs/F0);
%         bw = (F0/(Fs/2))/q;
%         [bComb, aComb] = iircomb(n, bw, 'notch');
%         bpRange = [1.5, 150];
%         [bBand, aBand] = butter(4, bpRange/(Fs/2));
%         Fs_cache = Fs;
%     end  
% 
%     emgProc = filtfilt(bComb, aComb, emgRaw);
%     emgProc = filtfilt(bBand, aBand, emgProc);
%     % emgDiff = emgProc(:,1) - emgProc(:,2);
% end


function trigger_labels = detect_triggers(triggerData,actualFs,subjId,sessionIdx,AV)

    % 找到每一个Trigchn的极差，从而判断有效Trigger通道

    % figure
    % 
    % plot((1:length(triggerData)),triggerData);
    trig_range = range(triggerData , 1);

    % 聚类分析极差显著较大的TrigChn
    k = 2; % 聚类数量

    % Use_trig 记录了有效的TrigChn
    [idx, ~] = kmeans(trig_range', k);

    if mean(trig_range(idx == 1)) >= mean(trig_range(idx == 2))

        Use_trig = triggerData(:,idx==1);

    else

        Use_trig = triggerData(:,idx==2);

    end



    % 下面寻找Trig的时间（先找Index，除以Fs得到时间）

    Difftrig = Use_trig(2:end,:)-Use_trig(1:end-1,:);

    Trig_cell = cell(size(Difftrig,2));


    % 对每一列的逐行差值进行聚类（根据数据特征大致可以分成4类，取显著大的数据作为Trigger_idx）

    for col = 1:size(Difftrig,2)

        Chn_to_analyze = Difftrig(:,col);

        k = 4;

        [idx, ~] = kmeans(Chn_to_analyze , k);

        % 找到均值最大的idx

        mean4 = zeros(1,4);

        for id = 1:4

            mean4(id) = mean(Chn_to_analyze(idx==id));

        end

        [~ , max_id] = max(mean4);

        Use_id = find(idx == max_id);

        Trig_cell{col} = Use_id;

        % Trig_cell{col} = find(Chn_to_analyze>max((Chn_to_analyze)+min(Chn_to_analyze))/2);
    end

    combined = [];

    for m = 1:size(Difftrig,2)

        combined = [combined; Trig_cell{m}];

    end

    sorted_combined = sort(combined, 'ascend');
    
    mean_diff = mean(sorted_combined(2:end)-sorted_combined(1:end-1));


    % 舍去相邻差值较小的异常值
    div = 2;     % 设定筛选阈值

    mean_diff = mean_diff / div;

    filtered = sorted_combined(1); 
    
    i = 1;
    % while true
    % 
    %     i = i+1;
    %     % 计算当前元素与前一个保留元素的差值
    %     if (sorted_combined(i) - filtered(end)) >= mean_diff  && (sorted_combined(i)-filtered(end)) < actualFs*11.5
    %         filtered = [filtered, sorted_combined(i)]; % 保留该元素
    % 
    %     elseif (sorted_combined(i)-filtered(end)) > 11.5*actualFs
    %         filtered = [filtered, round((sorted_combined(i)+filtered(end))/2)];
    %         i = i-1;
    % 
    %     end
    % 
    %     if i==length(sorted_combined)
    %         break;
    %     end
    % 
    % end


    while true

        i = i+1;
        % 计算当前元素与前一个保留元素的差值
        if (sorted_combined(i) - filtered(end)) >= mean_diff  
            filtered = [filtered, sorted_combined(i)]; % 保留该元素
        end

        % elseif (sorted_combined(i)-filtered(end)) > 11.5*actualFs
        %     filtered = [filtered, round((sorted_combined(i)+filtered(end))/2)];
        %     i = i-1;
        % 
        % end

        if i==length(sorted_combined)
            break;
        end

    end


   

    

    filtered_new = filtered(1);

   

    if subjId == 13 && sessionIdx == 1 && AV ==1 

        filtered_new = [];

        filtered_new(1,1) = 21154;

    elseif subjId == 13 && sessionIdx ==2 && AV ==1 

        filtered_new = [];

        filtered_new(1,1) = 24443;

    elseif subjId == 16 && sessionIdx == 1 && AV ==1 

        filtered_new = [];

        filtered_new(1,1) = 29375;

    elseif subjId == 16 && sessionIdx == 2 && AV ==1 

        filtered_new = [];

        filtered_new(1,1) = 14014;

    elseif subjId == 19 && sessionIdx ==1 && AV ==1 

        filtered_new = [];

        filtered_new(1,1) = 33564;

    elseif subjId == 19 && sessionIdx ==2 && AV ==1 

        filtered_new = [];

        filtered_new(1,1) = 17995;

    elseif subjId == 29 && sessionIdx ==1 && AV ==1 

        filtered_new = [];

        filtered_new(1,1) = 21641;

    elseif subjId == 29 && sessionIdx ==2 && AV ==1 

        filtered_new = [];

        filtered_new(1,1) = 26746;

    elseif subjId == 42 && sessionIdx ==2 && AV ==1

        filtered_new = [];

        filtered_new(1,1) = 34957;

    

    elseif subjId == 10 && sessionIdx == 1 && AV ==2

        filtered_new = [];

        filtered_new(1,1) = 1135;

    elseif subjId == 10 && sessionIdx ==2 && AV ==2

        filtered_new = [];

        filtered_new(1,1) = 16756; 

        

    elseif subjId == 14 && sessionIdx ==1 && AV ==2

        filtered_new = [];

        filtered_new(1,1) = 16214;

    elseif subjId == 14 && sessionIdx ==2 && AV ==2

        filtered_new = [];

        filtered_new(1,1) = 23612;
 
    elseif subjId == 16 && sessionIdx ==1 && AV ==2

        filtered_new = [];

        filtered_new(1,1) = 33709;

    elseif subjId == 16 && sessionIdx ==2 && AV ==2

        filtered_new = [];

        filtered_new(1,1) = 20021;

    elseif subjId == 18 && sessionIdx ==1 && AV ==2

        filtered_new = [];

        filtered_new(1,1) = 36946;

    elseif subjId == 18 && sessionIdx ==2 && AV ==2

        filtered_new = [];

        filtered_new(1,1) = 18222;

    end

    if AV == 2

        if subjId == 18 && sessionIdx == 1

            if filtered(end)+2*actualFs<435054

                filtered = [filtered,435054];
                
            end

        end

    elseif AV == 1

        if subjId == 42 && sessionIdx == 2

            if filtered(end)+2*actualFs<451528

                filtered = [filtered,451528];
                
            end

        end

    end

    % disp(filtered);


    for i = 1:length(filtered)

        if filtered(i) - filtered_new(end) > mean_diff && filtered(i) - filtered_new(end) < 12*actualFs

            filtered_new = [filtered_new;filtered(i)];

        elseif filtered(i) -filtered_new(end) > 16*actualFs && filtered(i) - filtered_new(end) < 20*actualFs

            filtered_new = [filtered_new;round((filtered_new(end)+filtered(i))/2);filtered(i)];

        elseif filtered(i) - filtered_new(end) >25*actualFs && filtered(i) - filtered_new(end) < 29*actualFs

            filtered_new = [filtered_new;round(filtered_new(end)*2/3+filtered(i)/3);round(filtered_new(end)/3+filtered(i)*2/3);filtered(i)];

        elseif filtered(i) - filtered_new(end) > 33*actualFs && filtered(i) - filtered_new(end) < 39*actualFs

            filtered_new = [filtered_new;round(filtered_new(end)*3/4+filtered(i)/4);round(filtered_new(end)/2+filtered(i)/2);round(filtered_new(end)/4+filtered(i)*3/4);filtered(i)];

        elseif filtered(i) - filtered_new(end) > 41*actualFs && filtered(i) - filtered_new(end) < 49*actualFs

            filtered_new = [filtered_new;round(filtered_new(end)*4/5+filtered(i)/5);round(filtered_new(end)*3/5+filtered(i)*2/5);round(filtered_new(end)*2/5+filtered(i)*3/5);round(filtered_new(end)*1/5+filtered(i)*4/5);filtered(i)];

        elseif filtered(i) - filtered_new(end) > 51*actualFs && filtered(i) - filtered_new(end) < 57*actualFs

            filtered_new = [filtered_new;round(filtered_new(end)*5/6+filtered(i)/6);round(filtered_new(end)*4/6+filtered(i)*2/6);round(filtered_new(end)/2+filtered(i)/2);round(filtered_new(end)*2/6+filtered(i)*4/6);round(filtered_new(end)/6+filtered(i)*5/6);filtered(i)];

        elseif filtered(i) - filtered_new(end) > 60*actualFs && filtered(i) - filtered_new(end) < 66*actualFs

            filtered_new = [filtered_new;round(filtered_new(end)*6/7+filtered(i)/7);round(filtered_new(end)*5/7+filtered(i)*2/7);round(filtered_new(end)*4/7+filtered(i)*3/7);round(filtered_new(end)*3/7+filtered(i)*4/7);round(filtered_new(end)*2/7+filtered(i)*5/7);round(filtered_new(end)/7+filtered(i)*6/7);filtered(i)];

        end

    end


   clear filtered

   filtered = filtered_new;


    disp(length(filtered));

    trigger_labels = zeros(size(triggerData,1),1);

    trigger_labels(filtered) = 1;

    

    % disp(length(filtered));

    % disp(filtered);



    
    % max_portion = 0.9;
    % 
    % 
    % 
    % threshold = max(triggerData,[],1)*max_portion + min(triggerData,[],1)*(1-max_portion);
    % 
    % midValues = threshold;
    % 
    % trigMask = (triggerData(2:end,:) >= midValues) & (triggerData(1:end-1,:) < midValues);
    % [rows, cols] = find(trigMask);
    % trigger_labels = zeros(size(triggerData,1),1);
    % trigger_labels(rows) = cols;
    % 
    % t = (1:length(triggerData))/500;
    % plot(t,triggerData);


    % midValues = (max(triggerData,[],1) + min(triggerData,[],1)) / 2;
    % disp(midValues);
    % trigMask = (triggerData(2:end,:) >= midValues) & (triggerData(1:end-1,:) < midValues);
    % [rows, cols] = find(trigMask);
    % trigger_labels = zeros(size(triggerData,1),1);
    % trigger_labels(rows) = cols;

end

function Data = cFilterD_EEG(data, Fs, k, OME)
    if k==1 || k==2
        Data = data;
        if length(OME) ~= 2
            error('The last variable shall be an 1X2 array, check the input!');
        end
        if k == 2
            % IIRCOMB filter       
            Fo = 50; q = 20;
            n = round(Fs/Fo);
            bw = (Fo/q)/(Fs/2); % 3dB bandwidth = Fo/q
            [B, A] = iircomb(n, bw, 'notch');
            Data = filtfilt(B, A, Data);
        end
        if Fs > OME(2)*2
            % Bandpass filter
            h = fdesign.bandpass('N,F3dB1,F3dB2', 8, OME(1), OME(2), Fs);
        else
            % Highpass filter
            h = fdesign.highpass('N,F3dB', 4, OME(1), Fs);
        end
        Hd = design(h, 'butter');
        Data = filtfilt(Hd.sosMatrix, Hd.scaleValues, Data);
    end
    if k==0
        if length(OME)>1
            error('This is low pass filter, OME length shall be 1!');
        end
    
        h  = fdesign.lowpass('N,F3dB', 4, OME, Fs);
        Hd = design(h, 'butter');
        Data = filtfilt(Hd.sosMatrix, Hd.scaleValues, data);
    end
end

function [Data, CHN] = CAR_reRef(config, SEEG, good_channels, subjId)
    EleCTX_dir = fullfile(config.EleCTX_dir, sprintf('P%d',subjId));
    % load(fullfile(EleCTX_dir, 'electrode_raw.mat'), 'elec_Info');
    load(fullfile(EleCTX_dir, 'SignalChanel_Electrode_Registration.mat'), 'CHN');
    if subjId == 38
        Sub_Chn = [64,64,64,30,20];

    elseif subjId == 43
        Sub_Chn = [64,64,32,24,24];
    end



    
    good_mask = false(size(CHN));
    good_mask(good_channels) = true;

    elec_number_end = cumsum(Sub_Chn);
    elec_number_begin = [1, elec_number_end(1:end-1) + 1];

    for jj = 1:length(Sub_Chn)

        elec_group = elec_number_begin(jj):elec_number_end(jj);
        Chn_Seq = CHN(elec_group);         % 组内通道        行：1      列：通道
        is_good = good_mask(Chn_Seq);      % 组内好通道判定     行：1     列：通道
        SEEG_shaft = SEEG(:, Chn_Seq);     % 组内选定通道数据   行：时间  列：通道
        SEEG_shaft_good = SEEG_shaft(:,is_good);
        to_minus = mean(SEEG_shaft_good,2);
        SEEG_shaft = SEEG_shaft - to_minus;
        SEEG(:,Chn_Seq) = SEEG_shaft;


    end

    Data = SEEG;










end





function [Data, CHN] = Laplacian_reRef(config, SEEG, good_channels, subjId)
    EleCTX_dir = fullfile(config.EleCTX_dir, sprintf('P%d', subjId));
    load(fullfile(EleCTX_dir, 'electrode_raw.mat'), 'elec_Info');
    load(fullfile(EleCTX_dir, 'SignalChanel_Electrode_Registration.mat'), 'CHN');
    
    Data = SEEG;
    good_mask = false(size(CHN));
    good_mask(good_channels) = true;

    Sub_Chn = cell2mat(elec_Info.number);
    elec_number_end = cumsum(Sub_Chn);
    elec_number_begin = [1, elec_number_end(1:end-1) + 1];
    
    for jj = 1:length(Sub_Chn)
        elec_group = elec_number_begin(jj):elec_number_end(jj);
        Chn_Seq = CHN(elec_group);
        is_good = good_mask(Chn_Seq);
        SEEG_shaft = SEEG(:, Chn_Seq);

        n = Sub_Chn(jj);
        coeff_left = zeros(1, n);
        coeff_right = zeros(1, n);
        

        if n >= 2 % should be more than 2 channels on 1 shaft
            if is_good(1) && is_good(2)
                coeff_right(1) = 1;
            end
            if is_good(end) && is_good(end-1)
                coeff_left(end) = 1;
            end
        end

        iVec = 1:n;
        mid_valid = (iVec > 1) & (iVec < n) & is_good; % boolean
        if any(mid_valid)
            valid_idx = find(mid_valid);
            left_valid = is_good(valid_idx - 1);
            right_valid = is_good(valid_idx + 1);

            both_valid = left_valid & right_valid;
            left_only = left_valid & ~both_valid;
            right_only = right_valid & ~both_valid;

            coeff_left(valid_idx) = 0.5*both_valid + 1*left_only;
            coeff_right(valid_idx) = 0.5*both_valid + 1*right_only;
        end

        SEEG_left = zeros(size(SEEG_shaft));        SEEG_right = SEEG_left;
        SEEG_left(:, 2:end) = SEEG_shaft(:, 1:end-1);
        SEEG_right(:, 1:end-1) = SEEG_shaft(:, 2:end);
        Data(:, Chn_Seq) = SEEG_shaft - coeff_left.*SEEG_left - coeff_right.*SEEG_right;
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













%% Test Area（测试结束记得注释）



% D = Datacell{1}(:,end);
% tr = find(D > 0);


% standard_1 = [39605,48803,58017,67248,76496,85777,95058,104306,113588,122885,132217,141514,150829,160143,169525,178972,188304,197618,206983,216330,225728,235143,244557,254005,263453,272934,282415,291880,301344,310809,320290,329754,339252,348833,358348,367879,377410,386908,396489,406053,415668,425699,435263,444844,454492]';
% standard_2 = [14809,23973,33154,42369,51584,60848,70129,79377,88658,97939,107237,116552,125966,135247,144579,153910,163241,172622,182087,191484,200832,210197,219595,229009,238424,247871,257303,266734,276182,285679,295277,304775,314256,323771,333318,342866,352364,361962,371509,381141,390722,400286,409884,419465,429080]';
% Z_1 = mean(abs(Trigger_ind_cell{1,1} - standard_1))/actualFs;
% Z_2 = mean(abs(Trigger_ind_cell{1,2} - standard_2))/actualFs;









% A = [1  13    20    25    11     8    17    22    29     5     9    28    19    16    18    10    15    17 ...
%  4     6    22    25    18    28    29    21    16    26    14     9    17    19    28    26    22     6  ...
%  16    11     1    27    12    10     2    14     9     2     6    24    26     3     4    10    21    12  ...
% 7    30    15    23    14    27    13    12    11     5     7     2     8    23     1     3     2    24 ...
% 19    27    20     3    13    23     5     7     4     8    18     9    25    20     9    21    15    24];