% 目前，运行出的gamma_epoch中包含部分Trigger（约占总数的5%）附近0.2s（由参数extent决定）的seeg数据决定，
% 第一维度为时间，第二维度为Trial，第三维度为Chn。
% 测试是否有显著响应的模块得到保留，可以更改平均值的倍数调整阈值设定。
% 已按照CHN顺序进行排序

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
extent1 = 1;     % 决定Trigger前多长的时间（单位：s）
extent2 = 8;     % 决定Trigger后多长的时间（单位：s）


for subjId = targetSubjects

    % con代表了第几个P

    count = 0; con = con+1;

    fprintf("Constructing gamma-epoch of P %d \n", subjId);

    % subInfo = get_subject_info(subjId);

    % Fs = subInfo.Fs;   actualFs = Fs;

    sessionNo = 1:2;

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
    
        g2 = size(datacell,2)-2;     % 有效通道数
    
        g3 = size(index_of_stimulus_onset,1);   % Trigger所代表的刺激数
    
        g4 = (extent1+extent2) * actualFs;               % 特定时间内的采样点数，4=2+2
    
        % gamma_epoch = zeros(g4 , g3 , g2);      % 维度1：g2->有效通道   g3->刺激试次   g4->刺激前后的一段时间
    
        for chn = 1:g2
    
            
    
            fprintf("Channel %d/%d is being processed.\n",chn,g2);
    
    
    
            selected_Chn = datacell(:,chn);
    
            seegDiff_smooth = smooth(abs(selected_Chn), 0.025*actualFs);
    
            for trial = 1:g3  % the segment number (i th)
    
                    count = count + 1;
    
                        % 开始构造Gamma3维矩阵
                        
                        
                        A = selected_Chn( ...
                            max(1,index_of_stimulus_onset(trial)-actualFs*extent1): ...
                            min(index_of_stimulus_onset(trial)+actualFs*extent2-1,length(Trig)));
        
                        
        
                        if length(A) == g4
        
                          gamma_epoch(:, trial, chn) = A;
        
   
                        end
                        

                    
    
            end
        end
    
        Gamma_epoch_cell{con,cons} = gamma_epoch;
    
        disp(size(gamma_epoch));


    end
end



