%% 确定哪些通道是我们希望保留的有响应的通道
% TF阶段必须使用全部通道





%% 输入：cwt_result_mean（cell）{Chn,gesture}（频率*时间）   f{Chn,gesture}（频率*时间）  subjId  
%% 输出：Chn_var



% 双尾p值检验样本量
number_of_repetitions = 100;



% con=0;
% 
% g1 = size(targetSubjects,2);
% 
% reactive_locations = cell(g1,2);
% 
% non_reactive_locations = cell(g1,2);
% 
% r_obs_cell = cell(g1,2);

Bandrange = cell(1,6);      % 5种频段的波
% Bandname = ['Delta','Theta','Alpha','Beta','Gamma','High-gamma'];
Bandrange{1,1} = [0.5,4];    % Delta
Bandrange{1,2} = [4,8];      % Theta
Bandrange{1,3} = [8,12];     % Alpha
Bandrange{1,4} = [12,30];    % Beta
Bandrange{1,5} = [30,80];    % Gamma
Bandrange{1,6} = [80,150];   % High-gamma

% Chn_react = cell(3,3); % Dim1：频带范围；Dim2：手势。

% check1 = 1*actualFs;   % 休息和想象的交界处
% 
% check2 = 4*actualFs;   % 想象和执行的交界处
% 
% check3 = 7*actualFs;   % 想象和休息的交界处

checkrange = cell(1,3);   % 3个检查点，检查点前后
% checkname = ['Rest & Imag','Rest and Exec'];
% 固定检验时间
checkrange{1,1} = 1:1*actualFs;      % 休息时间（0~1s）

checkrange{1,2} = (1*actualFs+1):3*actualFs;    % 运动想象时间（1~3s）

checkrange{1,3} = (3*actualFs+1):5*actualFs;    % 运动执行时间（4~6s）

Chn_num = size(cwt_result_mean,1);     % 通道数量

p_array = zeros(Chn_num,3,6,2);     % dim1:Chn_num  dim2:gesture  dim3:bandrange  dim4:checkpoint

% count_channels = zeros(1,Chn_num);

for chn = 1:Chn_num

    fprintf('Processing channel %d/%d \n',chn,Chn_num);

    can_be_adopted = 0;  % default: 该通道无法被纳入

    to_be_analyzed = cell(1,3);   

    for g = 1:3

        to_be_analyzed{1,g} = cwt_result_mean{chn,g};   % 1*3（gesture）的cell   每一个都是一个cwt_result的分析结果

    end
     
    % count = 0;     % 响应达标数量(对每个通道而言，总共有27处可检查，3*3*3）

    % 对每个手势单独分析，再对每个频带单独分析，如果p值检验p<0.05将该通道纳入备选通道

    for g = 1:3  % 3大手势

        fprintf("  Analyzing gesture %d \n",g);

        gesture_analysis = to_be_analyzed{1,g};    % cwt分析结果（矩阵）

        frequency_list = frequency_record{chn,g};

        for f = 1:6  % 6大频带

            fprintf("    Analyzing bandrange %d \n",f);

            frequency = frequency_list(frequency_list>Bandrange{1,f}(1) & frequency_list<Bandrange{1,f}(2));

            frequency_analysis = gesture_analysis(frequency_list>Bandrange{1,f}(1) & frequency_list<Bandrange{1,f}(2),:);

            for t = 1:2  % 2大检验点

               fprintf("      Analyzing checkpoint %d \n",t);

               checkpoint_analysis = mean(frequency_analysis);

               % before = checkpoint_analysis(checkrange{t,1});
               % 
               % after = checkpoint_analysis(checkrange{t,2});

               rest = checkpoint_analysis(checkrange{1,1});

               test = checkpoint_analysis(checkrange{1,t});


               N_rest = length(rest);

               N_test = length(test);

               % permutation テスト  （排列测试）

               

               reverseStr = '';

               labels = [zeros(1,N_rest),ones(1,N_test)];

               r_obs = corr([rest,test]',labels','Type','Spearman');

               r_pdf = zeros(1,number_of_repetitions);
               
               tmp = [rest,test];

               for repeat = 1:number_of_repetitions

                   w = tmp(randperm(length(tmp)));

                   r_pdf(repeat) = corr(w',labels','Type','Spearman');

               end

               % p_value = 1-normcdf(abs(r_obs), mean(r_pdf,2), std(r_pdf));

               % chiver = chi2gof(r_pdf);

               % 计算Z分数

               m_pdf = mean(r_pdf);

               s_pdf = std(r_pdf);

               z = (r_obs-m_pdf)/s_pdf;

               % 计算双尾P值

               p_value = 2*(1-normcdf(abs(z)));

               p_array(chn,g,f,t) = p_value;

            end
        end
    end
end

%% 

% p_array构成：dim1:chn   dim2:gesture(3)   dim3:frequency(6)
% dim4:t(2)



% usechn_idx = find(count_channels>20);

% usechn = Chn_wanted(usechn_idx);

% [~, name] = system('hostname');
% if strcmp(strip(name), 'ACE1999')
config.root_dir = fullfile('D:', 'Project','Code_Tutorial','03_Imaginary_Gesture_Coding','Coding_Gestures_Imaginary');
config.EleCTX_dir = fullfile(config.root_dir, '3_Brain_Electrodes');


EleCTX_dir = fullfile(config.EleCTX_dir, sprintf('P%d',subjId));

if subjId ~= 38 && subjId ~=43

    load(fullfile(EleCTX_dir, 'electrode_raw.mat'), 'elec_Info');

    Sub_Chn = cell2mat(elec_Info.number);

end



load(fullfile(EleCTX_dir,'electrodes_Final_Norm.mat'),'elec_Info_Final_wm');
Chn_info = elec_Info_Final_wm.name;       % Letters & numbers  % Maybe of no use due to the continuous property of elec_Info
name = elec_Info.name;    % Letters ONLY
name = cell2mat(name);
load(fullfile(EleCTX_dir, 'SignalChanel_Electrode_Registration.mat'), 'CHN');
if subjId == 38
    Sub_Chn = [64,64,64,30,20];
elseif subjId == 43
    Sub_Chn = [64,64,32,24,24];
end

Sub_Chn_add = cumsum(Sub_Chn);


% Now we get the argument Sub_Chn & CHN & Chn_info

% which means:
% Sub_Chn:  How many channels each part own (following the order of name)
% CHN:      Show the exact order
% Chn_info: Show the name of each channel

% Work to do next:
% 1. Sort the cahnnels into groups. [ ]
% 2. 

Alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'];

sign = 0;

starting_point = 0;

Chn_axis = zeros(1,Sub_Chn_add(end));    % 用于储存按字母顺序排列的通道序列（后续需要将坏通道去除）

BoundaryN = zeros(1:length(name));

BoundaryL = zeros(1:length(name));

for i = 1:length(Alphabet)

    sign = sign+1;

    Ind = find(name == Alphabet(sign));   % 对字母的指针

    while isempty(Ind)

        sign = sign+1;

        if sign>26 

            break;

        end

        Ind = find(name == Alpahebet(sign));

    end

    Letter = name(Ind);       % 指向字母

    if Ind~=1

        Chn_idx = CHN((Sub_Chn_add(Ind-1)+1):Sub_Chn_add(Ind));

    else

        Chn_idx = CHN(1:Sub_Chn_add(1));     % Chn_axis的一部分，记录了当下字母的通道

    end

    Chn_axis((starting_point+1):(starting_point+Sub_Chn(Ind))) = Chn_idx;

    BoundaryN(i) = starting_point+1;

    BoundaryL(i) = Letter;

    starting_point = starting_point + Sub_Chn(Ind);

end


%  Chn_axis包含坏通道

Chn_axis_good = [];

for i = 1:length(Chn_axis)

    if ismember(Chn_axis(i),good_channels)

        Chn_axis_good = [Chn_axis_good,Chn_axis(i)];

    else

        c_idx = find(BoundaryN>i);

        BoundaryN(c_idx) = BoundaryN(c_idx) -1;


    end


end


% Chn_axis_good保留了所有好通道并且按字母顺序排序。

% BoundaryN标记了所有不同字母交界处的序号

% BoundaryL是按顺序排列的字母列表

% p_array构成：dim1:chn   dim2:gesture(3)   dim3:frequency(6)
% dim4:t(2)


p_array_cell = cell(1,3);

for g = 1:3

    f = 1:150;

    [~,rank] = sort(Chn_axis_good);

    to_be_analyzed = p_array(:,g,:,:);

    to_be_analyzed = to_be_analyzed(rank,:,:);

    to_be_plotted = to_be_analyzed(:,:,2)./to_be_analyzed(:,:,1);    % 纵坐标为6个频段

    to_be_stored = zeros(150,length(Chn_axis_good));    % 纵坐标为150个频率值

    for i = 1:150

        if i>Bandrange{1,1}(1) && i<Bandrange{1,1}(2)

            to_be_stored(i,:) = to_be_plotted(:,1)';

        elseif i>Bandrange{1,2}(1) && i<Bandrange{1,2}(2)

            to_be_stored(i,:) = to_be_plotted(:,2)';

        elseif i>Bandrange{1,3}(1) && i<Bandrange{1,3}(2)

            to_be_stored(i,:) = to_be_plotted(:,3)';

        elseif i>Bandrange{1,4}(1) && i<Bandrange{1,4}(2)

            to_be_stored(i,:) = to_be_plotted(:,4)';

        elseif i>Bandrange{1,5}(1) && i<Bandrange{1,5}(2)

            to_be_stored(i,:) = to_be_plotted(:,5)';

        elseif i>Bandrange{1,6}(1) && i<Bandrange{1,6}(2)

            to_be_stored(i,:) = to_be_plotted(:,6)';

        end

    end

    % p_array_cell{1,g} = to_be_stored;

    % to_be_stored: 150*Chn_num（画图素材获得，还需卷积一下）

    core_size = 8;
                
    smooth_window = ones(core_size,1) / core_size;

    log_power = to_be_stored;

    smoothed_log_power = zeros(size(log_power));

    for chn = 1:size(log_power, 2)
        smoothed_log_power(:, chn) = conv(log_power(:,chn), smooth_window, 'same');
    end

    to_be_stored = smoothed_log_power;

    to_be_stored = 20*log10(to_be_stored);

    p_array_cell{1,g} = to_be_stored;


end


for g = 1:3

    figure;
    imagesc(1:length(Chn_axis_good),f,p_array_cell{1,g});
    hold on
    set(gca,"YDir","normal");
    % set(gca,'xtick',[],'ytick',[]);
    shading interp
    colormap("jet");
    axis tight;
    % axis xy;
    view(0, 90);
    xlabel('Channel');
    ylabel('Frequency (Hz)');
    ylim([0 150]);
    filename = sprintf('P%d Gesture%d',subjId,g);
    title(filename);
    colorbar;

    V = cwt_result_mean{Chn,gesture};

    % disp(prctile(V(:),95));
    clim([prctile(V(:),2) prctile(V(:),98)]);

    y_min = 0;

    y_max = max(f);

    for k = 1:length(BoundaryN)

        x_loc = BoundaryN(i);
        
        line([x_loc, x_loc], [y_min, y_max], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1);

        hold on


    end




end


