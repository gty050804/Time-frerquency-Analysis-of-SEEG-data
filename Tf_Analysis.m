% P17 129:139



% 由于目前在高频段上时间上的离散型较强，所以画图上尝试采用对数卷积平滑进行处理。
% 处理效果：提升

% 目前更正了原代码中出现的错误，再修正后得到类似理论图像，但可视化仍有待改进。


% 效果较好：P8   P9  P20  P12  
% 效果不好：P16  P13

% 电极通道映射
% Normalization   Okay
% 图片优化(colormap)，排版   
% (axes)

% 特征分类：低频成分   
% Location Identification：休息1s，想象2s；想象1s，执行2s。（挑选有响应的通道）
% Prediction


%% 输入：targetsubjexts 
%        Gamma_epoch
%        index_of_...
%        Trigger_ind_cell
%        g2
%        goodname

%% 输出：cwt_result_mean  （Chn，gesture）



% 目前仅支持单P分析


con=0;con1=0;
step = 1;    % 插值精度
% Channel to be analyzed
Chn_sel_1 = cell(length(targetSubjects),1);


% 1  现在还没改好
p = input("是否分析所有（若是，请输入1；若否，请输入0。）：");




if p==1

    Chn_wanted = 1:size(gamma_epoch_cell{1,1},3);

    Number_of_triggers = zeros(length(targetSubjects),3);   

    for con1 = 1:length(targetSubjects)
    
            Chn_sel_1{con1,1} = Chn_wanted;  

            
        
            % N(con1,1,gesture) = size(Gamma_epoch_cell{con1,1},2)+size(Gamma_epoch_cell{con1,2},2);

            Number_of_triggers(con1,1) = length(find(index_of_stimulus_onset_cell{con1,1}<=Trigger_ind_cell{con1,1}(15)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}<=Trigger_ind_cell{con1,2}(15)));

            Number_of_triggers(con1,2) = length(find(index_of_stimulus_onset_cell{con1,1}>Trigger_ind_cell{con1,1}(15)&index_of_stimulus_onset_cell{con1,1}<=Trigger_ind_cell{con1,1}(30)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}>Trigger_ind_cell{con1,2}(15)&index_of_stimulus_onset_cell{con1,2}<=Trigger_ind_cell{con1,2}(30)));

            Number_of_triggers(con1,3) = length(find(index_of_stimulus_onset_cell{con1,1}>Trigger_ind_cell{con1,1}(30)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}>Trigger_ind_cell{con1,2}(30)));




    end
end



if p~=1

    fprintf("通道数量上限：%d \n",g2);

    Chn_wanted = input("请输入待分析通道：");

    Number_of_triggers = zeros(length(targetSubjects),3);   

    for con1 = 1:length(targetSubjects)
    
            Chn_sel_1{con1,1} = Chn_wanted;  
      
            % N(con1,1,gesture) = size(Gamma_epoch_cell{con1,1},2)+size(Gamma_epoch_cell{con1,2},2);

            Number_of_triggers(con1,1) = length(find(index_of_stimulus_onset_cell{con1,1}<=Trigger_ind_cell{con1,1}(15)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}<=Trigger_ind_cell{con1,2}(15)));

            Number_of_triggers(con1,2) = length(find(index_of_stimulus_onset_cell{con1,1}>Trigger_ind_cell{con1,1}(15)&index_of_stimulus_onset_cell{con1,1}<=Trigger_ind_cell{con1,1}(30)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}>Trigger_ind_cell{con1,2}(15)&index_of_stimulus_onset_cell{con1,2}<=Trigger_ind_cell{con1,2}(30)));

            Number_of_triggers(con1,3) = length(find(index_of_stimulus_onset_cell{con1,1}>Trigger_ind_cell{con1,1}(30)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}>Trigger_ind_cell{con1,2}(30)));


    end
end
    
    

outputFolder = sprintf('D:\\Project\\Code_Tutorial\\03_Imaginary_Gesture_Coding\\Coding_Gestures_Imaginary\\CWT_data\\P%d',subjId);
out = cell(1,3);
out{1,1} = fullfile(outputFolder,'gesture1');
out{1,2} = fullfile(outputFolder,'gesture2');
out{1,3} = fullfile(outputFolder,'gesture3');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end 
if ~exist(out{1,1}, 'dir')
    mkdir(out{1,1});
end 
if ~exist(out{1,2}, 'dir')
    mkdir(out{1,2});
end 
if ~exist(out{1,3}, 'dir')
    mkdir(out{1,3});
end 

whole_progress = waitbar(0,"Don't worry. It's almost there!",'Name','TF_Analysis');
progress = 0;
Number_whole = 0;
for i = 1:size(Gamma_epoch_cell,1)
    for j = 1:size(Gamma_epoch_cell,2)
        Number_whole = Number_whole + size(Gamma_epoch_cell{i,j},2);
    end
end
Number_whole = Number_whole * length(Chn_wanted);


% window_size = 18;   % 选取18s的时间段作为时频分析的素材
% rng(42);      % 随机数种子

%% Time-Frequency Analysis
for subjId = targetSubjects

    clear cwt_result_mean
    % cwt_result_mean = cell(length(Chn_sel),1);

    % for u = 1:length(Chn_sel)
    % 
    %     cwt_result_mean{u,1} = [];
    % 
    % end
    % 
  

    subInfo = get_subject_info(subjId);

    % Fs = subInfo.Fs;   actualFs = Fs;

    con = con+1;

    cwt_result_mean = cell(length(Chn_sel_1{con,1}),3);


    % N = size(Datacell{1},1);
    % 
    % t = (0:N-1)/Fs;

    sessionNo = subInfo.Session_num;

    for cons = sessionNo

        fprintf("Analysing time-frequency of P%d Session%d... \n", subjId,cons);

        index_of_stimulus_onset = index_of_stimulus_onset_cell{con,cons};

        Chn_sel = Chn_sel_1{con,1};

        gamma_interest = Gamma_epoch_cell{con,cons};

        Chn_interest = gamma_interest(:,:,Chn_sel);

        m = 0;

        % fprintf('该试次下总Trigger数为%d \n',size(Chn_interest,2));
    
        % Trig_idx = input('请输入待分析的Trigger序号：');

        Trig_idx = 1:size(Chn_interest,2);

        
        % N = N+size(Chn_interest,2);

        for Chn = Chn_sel

            fprintf('Analysing Channel %d \n',Chn);

            % barname = sprintf('P%d Session%d Channel %d \n',subjId,cons,Chn);
            % 
            m = m+1;
            % 
            % Chn_progress = sprintf('Channel %d/%d',m,length(Chn_sel));
            % 
            % h = waitbar(0, Chn_progress, 'Name', barname);


            % to_be_analyzed = Chn_interest(:,Trig_idx,m);

            % t = (1:window_size*actualFs)/actualFs;
            % 
            % start_point = randi([1,(length(to_be_analyzed)-window_size*actualFs)]);
            % 
            % end_point = start_point+window_size*actualFs;
            % 
            % t = t()

            for i = 1:length(Trig_idx)

                to_be_analyzed = Chn_interest(:,Trig_idx,m);

                to_be_analyzed = to_be_analyzed(:,i);

                % to_be_analyzed = [to_be_analyzed(2*actualFs+1:end);to_be_analyzed(1:2*actualFs)];

                % 现在的to_be_analyzed已经是一个9s的行向量了。
                % 长度：9*actualFs
                % 下面进行点对称延拓（防止小波变换的边界效应）

                comp = flip(to_be_analyzed);

                comp1 = 2*to_be_analyzed(1) - comp;

                comp2 = 2*to_be_analyzed(end) - comp;

                to_be_analyzed = [comp1;to_be_analyzed;comp2];

                t = (1:length(to_be_analyzed))/actualFs;

                [cwt_result,f] = cwt(to_be_analyzed,'morse',actualFs);

               % frequency_record{m,}

       

                cwt_result = abs(cwt_result);

                start_point = length(comp1)+1;

                end_point = length(comp1)+length(comp);

                cwt_result = cwt_result(:,start_point:end_point);

                 % 插值处理

                f_result = (1:(150/step))*step;

                f_result = flip(f_result);

                cwt_new = zeros(150/step,9*actualFs);    % 统一化后的cwt

                for col = 1:size(cwt_result, 2)
    
                     cwt_new(:, col) = interp1(f, cwt_result(:, col), f_result, 'linear', 'extrap');

                end

                cwt_result = cwt_new;

                norm_standard = [cwt_result(:,0.5*actualFs+1:1*actualFs),cwt_result(:,8*actualFs+1:9*actualFs)];

                cwt_result = cwt_result./repmat(median(cwt_result,2) ...
                    ,1,size(cwt_result,2));

                

                if index_of_stimulus_onset(i) <= Trigger_ind_cell{con,cons}(15)

                   

                    if isempty(cwt_result_mean{m,1})

                    cwt_result_mean{m,1} = cwt_result;

                    else
                    
                    cwt_result_mean{m,1} = cwt_result_mean{m,1}+cwt_result;

                    end

                    



                elseif index_of_stimulus_onset(i) > Trigger_ind_cell{con,cons}(15)  && index_of_stimulus_onset(i) <= Trigger_ind_cell{con,cons}(30)
                    
                   
                    if isempty(cwt_result_mean{m,2})

                    cwt_result_mean{m,2} = cwt_result;

                    else
                    
                    cwt_result_mean{m,2} = cwt_result_mean{m,2}+cwt_result;

                    end

                    

                else
                    
                    

                    if isempty(cwt_result_mean{m,3})

                    cwt_result_mean{m,3} = cwt_result;

                    else
                    
                    cwt_result_mean{m,3} = cwt_result_mean{m,3}+cwt_result;

                    end

                    

                end
                progress = progress+1;
                

                % waitbar(i/length(Trig_idx),h);
                waitbar(progress/Number_whole,whole_progress);

            end

            % close(h);
            

        end



    end

for gesture = 1:3

    for Chn = 1:length(Chn_sel)

                cwt_result_mean{Chn,gesture} = cwt_result_mean{Chn,gesture}/Number_of_triggers(con1,gesture);

                % 创建高斯滤波器
                sigma = 10;          % 高斯核的标准差，可以调整
                filterSize = [3 3];  % 滤波器大小
                
                % 应用高斯滤波
                cwt_result_mean{Chn,gesture} = imgaussfilt(cwt_result_mean{Chn,gesture}, sigma, 'FilterSize', filterSize);
                % cwt_result_mean{Chn,gesture} = imgaussfilt(cwt_result_mean{Chn,gesture}, sigma, 'FilterSize', filterSize);
                
                cwt_result_mean{Chn,gesture} = 20 * log10(cwt_result_mean{Chn,gesture});

                % Log-smoothing
                
                core_size = 500;
                
                smooth_window = ones(1, core_size) / core_size;

                log_power = cwt_result_mean{Chn,gesture};

                smoothed_log_power = zeros(size(log_power));

                for freq_idx = 1:size(log_power, 1)
                    smoothed_log_power(freq_idx, :) = conv(log_power(freq_idx, :), smooth_window, 'same');
                end

                cwt_result_mean{Chn,gesture} = smoothed_log_power;











                set(0,'DefaultFigureVisible', 'on');
                t = t(1:end_point-start_point+1);

                figure;
                imagesc(t,f_result,cwt_result_mean{Chn,gesture});
                hold on
                set(gca,"YDir","normal");
                shading interp
                colormap("jet");
                axis tight;
                % axis xy;
                view(0, 90);
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                ylim([0 150]);
                filename = sprintf('P%d Channel %s Gesture%d',subjId,goodname{Chn_sel(Chn)},gesture);
                title(filename);
                colorbar;
                % clim([0 prctile(cwt_result,95,'all')]);




                V = cwt_result_mean{Chn,gesture};

                % disp(prctile(V(:),95));
                clim([prctile(V(:),10) prctile(V(:),98)]);

                x_loc1 = round(size(cwt_result_mean{Chn,gesture},2)/9)/actualFs;

                x_loc2 = round(size(cwt_result_mean{Chn,gesture},2)/9*4)/actualFs;

                x_loc3 = round(size(cwt_result_mean{Chn,gesture},2)/9*7)/actualFs; 

                y_min = min(f);

                y_max = max(f);

                line([x_loc1, x_loc1], [y_min, y_max], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1);

                hold on

                line([x_loc2, x_loc2], [y_min, y_max], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1);

                hold on

                line([x_loc3, x_loc3], [y_min, y_max], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1);



                first_character = goodname{Chn_sel(Chn)}(1);

                outputFolder = fullfile(out{1,gesture},first_character);

                if ~exist(outputFolder, 'dir')
                    % 创建文件夹
                    mkdir(outputFolder);
                    
                end


                
                if Chn == length(Chn_sel)
                    fprintf("Saving figure of Channel%d. \n",Chn);
                else
                    fprintf("Saving figure of Channel%d... \n",Chn);
                end

                saveas(gcf, fullfile(outputFolder, filename),'png');





    end


end

end
close(whole_progress);



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



