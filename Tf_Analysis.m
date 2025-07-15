% 效果较好：P11 1;P11 2;




con = 0;con1=0;
% Channel to be analyzed
Chn_sel_1 = cell(length(targetSubjects),1);

p = input("是否分析所有（若是，请输入1；若否，请输入0。）：");


if p==1

    for con1 = 1:length(targetSubjects)

        Chn_sel_1{con1,1} = 1:size(Gamma_epoch_cell{con1,1},3);

    end

end



if p~=1

    Chn_wanted = input("请输入待分析通道：");

    for con1 = 1:length(targetSubjects)
    
            Chn_sel_1{con1,1} = Chn_wanted;
    
    end

end






% window_size = 18;   % 选取18s的时间段作为时频分析的素材
% rng(42);      % 随机数种子


%% Time-Frequency Analysis
for subjId = targetSubjects



    

    subInfo = get_subject_info(subjId);

    Fs = subInfo.Fs;   actualFs = Fs;

    con = con+1;

    % N = size(Datacell{1},1);
    % 
    % t = (0:N-1)/Fs;

    sessionNo = subInfo.Session_num;

    for cons = sessionNo

        fprintf("Analysing time-frequency of P%d Session%d... \n", subjId,cons);

        Chn_sel = Chn_sel_1{con,1};

        gamma_interest = Gamma_epoch_cell{con,cons};

        Chn_interest = gamma_interest(:,:,Chn_sel);

        m = 0;

        fprintf('该试次下总Trigger数为%d \n',size(Chn_interest,2));
    
        Trig_idx = input('请输入待分析的Trigger序号：');

        for Chn = Chn_sel

            

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

                % 现在的to_be_analyzed已经是一个9s的行向量了。
                % 长度：9*actualFs

                % 下面进行点对称延拓（防止小波变换的边界效应）

                comp = flip(to_be_analyzed);

                comp1 = to_be_analyzed(1) - comp;

                comp2 = to_be_analyzed(end) - comp;

                to_be_analyzed = [comp1;to_be_analyzed;comp2];

                t = (1:length(to_be_analyzed))/actualFs;

               [cwt_result,f] = cwt(to_be_analyzed,'morse',actualFs);

                cwt_result = abs(cwt_result);

                start_point = length(comp1)+1;

                end_point = length(comp1)+length(comp);

                cwt_result = cwt_result(:,start_point:end_point);

                cwt_result = cwt_result./repmat(median(cwt_result(:,1:end_point-start_point+1),2) ...
                    ,1,size(cwt_result,2));

                % disp(size(repmat(median(cwt_result(:,1:end_point-start_point+1),2) ...
                %     ,1,size(cwt_result,2))));

                t = t(1:end_point-start_point+1);

                % 绘制小波变换的时频图
                
                set(0,'DefaultFigureVisible', 'on');
                
                figure;
                surface(t, f, cwt_result, 'EdgeColor', 'none');
                axis tight;
                view(0, 90);
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                ylim([0 150]);
                filename = sprintf('P%d Session%d Channel%d Trigger%d',subjId,cons,Chn,Trig_idx(i));
                title(filename);
                colorbar;
                clim([0 prctile(cwt_result,99,'all')]);

               

                    outputFolder = 'D:\Project\Code_Tutorial\03_Imaginary_Gesture_Coding\Coding_Gestures_Imaginary\1_Raw_Data_All\Output\TF_Analysis\';
                    if Chn == Chn_sel(end)
                        fprintf("Saving figure of Channel%d. \n",Chn);
                    else
                        fprintf("Saving figure of Channel%d... \n",Chn);
                    end
    
                    % saveas(gcf, fullfile(outputFolder, filename),'png');

                





            end

            




        end


        










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



