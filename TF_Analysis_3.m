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

%% 输出：cwt_result_mean  （Chn，gesture）



% 目前仅支持单P分析


con=0;con1=0;
% Channel to be analyzed
Chn_sel_1 = cell(length(targetSubjects),1);


% 1  现在还没改好
p = input("是否分析所有（若是，请输入1；若否，请输入0。）：");




if p==1

    Chn_wanted = 1:size(gamma_epoch_cell{1,1},3);

    Number_of_triggers = zeros(length(targetSubjects),1,3);   

    for con1 = 1:length(targetSubjects)
    
            Chn_sel_1{con1,1} = Chn_wanted;  

            
        
            % N(con1,1,gesture) = size(Gamma_epoch_cell{con1,1},2)+size(Gamma_epoch_cell{con1,2},2);

            Number_of_triggers(con1,1,1) = length(find(index_of_stimulus_onset_cell{con1,1}<=Trigger_ind_cell{con1,1}(15)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}<=Trigger_ind_cell{con1,2}(15)));

            Number_of_triggers(con1,1,2) = length(find(index_of_stimulus_onset_cell{con1,1}>Trigger_ind_cell{con1,1}(15)&index_of_stimulus_onset_cell{con1,1}<=Trigger_ind_cell{con1,1}(30)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}>Trigger_ind_cell{con1,2}(15)&index_of_stimulus_onset_cell{con1,2}<=Trigger_ind_cell{con1,2}(30)));

            Number_of_triggers(con1,1,3) = length(find(index_of_stimulus_onset_cell{con1,1}>Trigger_ind_cell{con1,1}(30)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}>Trigger_ind_cell{con1,2}(30)));




    end
end



if p~=1

    Chn_wanted = input("请输入待分析通道：");

    Number_of_triggers = zeros(length(targetSubjects),1,3);   

    for con1 = 1:length(targetSubjects)
    
            Chn_sel_1{con1,1} = Chn_wanted;  

            
        
            % N(con1,1,gesture) = size(Gamma_epoch_cell{con1,1},2)+size(Gamma_epoch_cell{con1,2},2);

            Number_of_triggers(con1,1,1) = length(find(index_of_stimulus_onset_cell{con1,1}<=Trigger_ind_cell{con1,1}(15)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}<=Trigger_ind_cell{con1,2}(15)));

            Number_of_triggers(con1,1,2) = length(find(index_of_stimulus_onset_cell{con1,1}>Trigger_ind_cell{con1,1}(15)&index_of_stimulus_onset_cell{con1,1}<=Trigger_ind_cell{con1,1}(30)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}>Trigger_ind_cell{con1,2}(15)&index_of_stimulus_onset_cell{con1,2}<=Trigger_ind_cell{con1,2}(30)));

            Number_of_triggers(con1,1,3) = length(find(index_of_stimulus_onset_cell{con1,1}>Trigger_ind_cell{con1,1}(30)))...
            +length(find(index_of_stimulus_onset_cell{con1,2}>Trigger_ind_cell{con1,2}(30)));




    end
end
    
    

    






% window_size = 18;   % 选取18s的时间段作为时频分析的素材
% rng(42);      % 随机数种子


%% Time-Frequency Analysis
for subjId = targetSubjects


    % cwt_result_mean = cell(length(Chn_sel),1);

    % for u = 1:length(Chn_sel)
    % 
    %     cwt_result_mean{u,1} = [];
    % 
    % end
    % 
    % 

   
    

    subInfo = get_subject_info(subjId);

    Fs = subInfo.Fs;   actualFs = Fs;

    con = con+1;

    cwt_result_mean = cell(length(Chn_sel_1{con,1}),3);

    frequency_record = cell(length(Chn_sel_1{con,1}),3);

    

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

            m = m+1;
            
            

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

                to_be_analyzed = [to_be_analyzed(2*actualFs+1:end);to_be_analyzed(1:2*actualFs)];

                % 现在的to_be_analyzed已经是一个9s的行向量了。
                % 长度：9*actualFs
                % 下面进行点对称延拓（防止小波变换的边界效应）

                comp = flip(to_be_analyzed);

                comp1 = to_be_analyzed(1) - comp;

                comp2 = to_be_analyzed(end) - comp;

                to_be_analyzed = [comp1;to_be_analyzed;comp2];

                t = (1:length(to_be_analyzed))/actualFs;

               [cwt_result,f] = cwt(to_be_analyzed,'morse',actualFs);

               % frequency_record{m,}






                cwt_result = abs(cwt_result);

                start_point = length(comp1)+1;

                end_point = length(comp1)+length(comp);

                cwt_result = cwt_result(:,start_point:end_point);

                cwt_result = cwt_result./repmat(median(cwt_result(:,1*actualFs:3*actualFs),2) ...
                    ,1,size(cwt_result,2));

                

                if index_of_stimulus_onset(i) <= Trigger_ind_cell{con,cons}(15)

                    if isempty(cwt_result_mean{m,1})

                    cwt_result_mean{m,1} = cwt_result;

                    else
                    
                    cwt_result_mean{m,1} = cwt_result_mean{m,1}+cwt_result;

                    end

                    if isempty(frequency_record{m,1})

                    frequency_record{m,1} = f;

                    else

                    frequency_record{m,1} = frequency_record{m,1}+f;

                    end



                elseif index_of_stimulus_onset(i) > Trigger_ind_cell{con,cons}(15)  && index_of_stimulus_onset(i) <= Trigger_ind_cell{con,cons}(30)

                    if isempty(cwt_result_mean{m,2})

                    cwt_result_mean{m,2} = cwt_result;

                    else
                    
                    cwt_result_mean{m,2} = cwt_result_mean{m,2}+cwt_result;

                    end

                    if isempty(frequency_record{m,2})

                    frequency_record{m,2} = f;

                    else

                    frequency_record{m,2} = frequency_record{m,2}+f;

                    end

                else

                    if isempty(cwt_result_mean{m,3})

                    cwt_result_mean{m,3} = cwt_result;

                    else
                    
                    cwt_result_mean{m,3} = cwt_result_mean{m,3}+cwt_result;

                    end

                    if isempty(frequency_record{m,3})

                    frequency_record{m,3} = f;

                    else

                    frequency_record{m,3} = frequency_record{m,3}+f;

                    end

                end



                % if isempty(cwt_result_mean{m,1})
                % 
                %     cwt_result_mean{m,1} = cwt_result;
                % 
                % else
                % 
                %     cwt_result_mean{m,1} = cwt_result_mean{m,1}+cwt_result;
                % 
                % end



                % if i ~= 1
                % 
                %     cwt_result_mean = cwt_result_mean +cwt_result;
                % 
                % end
                % 
                % if i == 1
                % 
                %     cwt_result_mean = cwt_result;
                % 
                % end

                % disp(size(repmat(median(cwt_result(:,1:end_point-start_point+1),2) ...
                %     ,1,size(cwt_result,2))));

                % t = t(1:end_point-start_point+1);

                % 绘制小波变换的时频图
                
                % set(0,'DefaultFigureVisible', 'on');
                % 
                % 
                % figure;
                % surface(t, f, cwt_result, 'EdgeColor', 'none');
                % axis tight;
                % view(0, 90);
                % xlabel('Time (s)');
                % ylabel('Frequency (Hz)');
                % ylim([0 150]);
                % filename = sprintf('P%d Session%d Channel%d Trigger%d',subjId,cons,Chn,Trig_idx(i));
                % title(filename);
                % colorbar;
                % clim([0 prctile(cwt_result,99,'all')]);
                % 
                % 
                % 
                %     outputFolder = 'D:\Project\Code_Tutorial\03_Imaginary_Gesture_Coding\Coding_Gestures_Imaginary\1_Raw_Data_All\Output\TF_Analysis\';
                %     if Chn == Chn_sel(end)
                %         fprintf("Saving figure of Channel%d. \n",Chn);
                %     else
                %         fprintf("Saving figure of Channel%d... \n",Chn);
                %     end
                % 
                %     % saveas(gcf, fullfile(outputFolder, filename),'png');
                % 
                % 
                % 




            end

             % cwt_result_mean = cwt_result_mean/N(con,1);
             % 
             %    set(0,'DefaultFigureVisible', 'on');
             %    t = t(1:end_point-start_point+1);
             % 
             %    figure;
             %    surface(t, f, cwt_result_mean, 'EdgeColor', 'none');
             %    axis tight;
             %    view(0, 90);
             %    xlabel('Time (s)');
             %    ylabel('Frequency (Hz)');
             %    ylim([0 150]);
             %    filename = sprintf('P%d Channel%d',subjId,Chn);
             %    title(filename);
             %    colorbar;
             %    clim([0 prctile(cwt_result,99,'all')]);
             % 
             % 
             % 
             %        outputFolder = 'D:\Project\Code_Tutorial\03_Imaginary_Gesture_Coding\Coding_Gestures_Imaginary\1_Raw_Data_All\Output\TF_Analysis\';
             %        if Chn == Chn_sel(end)
             %            fprintf("Saving figure of Channel%d. \n",Chn);
             %        else
             %            fprintf("Saving figure of Channel%d... \n",Chn);
             %        end
             % 




        end



    end

for gesture = 1:3

    for Chn = 1:length(Chn_sel)

                cwt_result_mean{Chn,gesture} = cwt_result_mean{Chn,gesture}/Number_of_triggers(con1,1,gesture);

                % 创建高斯滤波器
                sigma = 10;          % 高斯核的标准差，可以调整
                filterSize = [3 3];  % 滤波器大小
                
                % 应用高斯滤波
                cwt_result_mean{Chn,gesture} = imgaussfilt(cwt_result_mean{Chn,gesture}, sigma, 'FilterSize', filterSize);
                cwt_result_mean{Chn,gesture} = imgaussfilt(cwt_result_mean{Chn,gesture}, sigma, 'FilterSize', filterSize);
                
                cwt_result_mean{Chn,gesture} = log10(cwt_result_mean{Chn,gesture});

                set(0,'DefaultFigureVisible', 'on');
                t = t(1:end_point-start_point+1);

                figure;
                imagesc(t, f, cwt_result_mean{Chn,gesture});
                set(gca,"YDir","normal");
                shading interp
                colormap("jet");
                axis tight;
                view(0, 90);
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                ylim([0 150]);
                filename = sprintf('P%d Channel %d Gesture%d',subjId,Chn_sel(Chn),gesture);
                title(filename);
                colorbar;
                % clim([0 prctile(cwt_result,95,'all')]);
                V = cwt_result_mean{Chn,gesture};
                % disp(prctile(V(:),95));
                clim([prctile(V(:),8) prctile(V(:),98)]);

                
                



                outputFolder = sprintf('D:\\Project\\Code_Tutorial\\03_Imaginary_Gesture_Coding\\Coding_Gestures_Imaginary\\1_Raw_Data_All\\Output\\TF_Analysis\\P%d',subjId);
                if ~exist(outputFolder, 'dir')
                    mkdir(outputFolder);
                end
                if Chn == Chn_sel(end)
                    fprintf("Saving figure of Channel%d. \n",Chn);
                else
                    fprintf("Saving figure of Channel%d... \n",Chn);
                end

                % saveas(gcf, fullfile(outputFolder, filename),'png');


                % figure;
                % for number = 1:4
                %     plot(t,cwt_result_mean{Chn,1}(number*5,:));
                %     hold on
                % end
                % title('Magnitude change with t');

             




    



    end







% 
% [cwt_result,f] = cwt(to_be_analyzed,actualFs);
% 
%                 cwt_result = abs(cwt_result);
% 
%                 start_point = length(comp1)+1;
% 
%                 end_point = length(comp1)+length(comp);

% 
% [cwt_result,f] = timefreq(to_be_analyzed,actualFs,'freqs',[0,150],'tlimits',[0 9*Fs],'winsize',500);
% 
%                 cwt_result = abs(cwt_result);
% 
%                 start_point = round(size(cwt_result,2)/3)+1;
% 
%                 end_point = round(size(cwt_result,2)/3*2)+1;







    % % 绘制小波变换的时频图
    % figure;
    % surface(t, f, abs(cwt_SEEG), 'EdgeColor', 'none');
    % axis tight;
    % view(0, 90);
    % xlabel('Time (s)');
    % ylabel('Frequency (Hz)');
    % title('CWT Magnitude');
    % colorbar;

    



    % 生成一个非平稳信号
    % fs = 1000;          % 采样频率 (Hz)
    % t = 0:1/fs:2;       % 时间向量
    % x = cos(2*pi*50*t) + cos(2*pi*120*t).*(t > 0.5 & t < 1.5); % 混合信号

    % 使用小波变换计算信号的时频表示
    % [cwt_result, f] = cwt(x, fs);
    % 
    % % 绘制小波变换的时频图
    % figure;
    % surface(t, f, abs(cwt_result), 'EdgeColor', 'none');
    % axis tight;
    % view(0, 90);
    % xlabel('Time (s)');
    % ylabel('Frequency (Hz)');
    % title('CWT Magnitude');
    % colorbar;
    % 

                % cwt_result_mean = cwt_result_mean/N(con,1);
                % 
                % set(0,'DefaultFigureVisible', 'on');
                % t = t(1:end_point-start_point+1);
                % 
                % figure;
                % surface(t, f, cwt_result_mean, 'EdgeColor', 'none');
                % axis tight;
                % view(0, 90);
                % xlabel('Time (s)');
                % ylabel('Frequency (Hz)');
                % ylim([0 150]);
                % filename = sprintf('P%d Channel%d',subjId,Chn);
                % title(filename);
                % colorbar;
                % clim([0 prctile(cwt_result,99,'all')]);
                % 
                % 
                % 
                %     outputFolder = 'D:\Project\Code_Tutorial\03_Imaginary_Gesture_Coding\Coding_Gestures_Imaginary\1_Raw_Data_All\Output\TF_Analysis\';
                %     if Chn == Chn_sel(end)
                %         fprintf("Saving figure of Channel%d. \n",Chn);
                %     else
                %         fprintf("Saving figure of Channel%d... \n",Chn);
                %     end
                % 
                    % saveas(gcf, fullfile(outputFolder, filename),'png');













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



