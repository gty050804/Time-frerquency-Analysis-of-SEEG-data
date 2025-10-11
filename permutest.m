% -[x] 挑选通道进行测试
% 在测试中发现计算开销大，且带通滤波后的信号似乎未能正确代表特定频段的信号强度

% -[ ] for循环过多，抛弃cell结构进行并行运算


% 需要用到的变量
% 1. Gamma_epoch_cell
% 2. good_channels

Chn_want = input("请输入待分析通道：");

use_channels = good_channels(Chn_want);




Bandrange = cell(1,6);      % 5种频段的波
% Bandname = ['Delta','Theta','Alpha','Beta','Gamma','High-gamma'];
Bandrange{1,1} = [0.5,4];    % Delta
Bandrange{1,2} = [4,8];      % Theta
Bandrange{1,3} = [8,12];     % Alpha
Bandrange{1,4} = [12,30];    % Beta
Bandrange{1,5} = [30,80];    % Gamma
Bandrange{1,6} = [80,150];   % High-gamma

number_of_repetitions = 50;
% 构建测试材料矩阵
sozai = zeros(length(use_channels),size(Gamma_epoch_cell{1,1},2)+size(Gamma_epoch_cell{1,2},2), ...
    length(Bandrange),3);    % 4-dim array   dim1:Chn  dim2:trigger  dim3:band  dim4:3
Gamma_cell = cell(1,sessionNum);    % 已挑选好通道，goodname与之一一对应
for i = 1:sessionNum
    Gamma_cell{1,i} = Gamma_epoch_cell{1,i}(:,:,use_channels);
end




for i = 1:sessionNum
    t_num = size(Gamma_cell{1,i},2);     % trigger数量
    ana_array = Gamma_cell{1,i};
    offset = mod(i-1,2)*size(Gamma_cell{1,1},2);
    
    % 预分配临时存储数组
    temp_sozai = zeros(length(use_channels), t_num, length(Bandrange), 3);
    
    parfor chn = 1:length(use_channels)
        chn_temp_sozai = zeros(t_num, length(Bandrange), 3);
        
        for trigger = 1:t_num
            to_be_analyzed = ana_array(:,trigger,chn);     % 待分析信号
            fprintf("Analysing trigger %d/%d \n",trigger,t_num);
            
            % 为每个band预分配结果
            band_results = zeros(length(Bandrange), 3);
            
            for band = 1:length(Bandrange)
                analyzed = bandpass(to_be_analyzed,[Bandrange{1,band}(1),Bandrange{1,band}(2)],actualFs);  % 带通滤波后的信号

                m1 = mean(analyzed(1:1*actualFs));    % 休息段均值
                m2 = mean(analyzed((1*actualFs+1):3*actualFs));   % 想象段均值
                m3 = mean(analyzed((4*actualFs+1):6*actualFs));   % 执行段均值
                
                band_results(band, :) = [m1, m2, m3];
            end
            
            chn_temp_sozai(trigger, :, :) = band_results;
        end
        
        % 将当前通道的结果存入临时数组
        temp_sozai(chn, :, :, :) = chn_temp_sozai;
    end
    
    % 将临时数组复制到最终结果中
    sozai(:, (offset+1):(offset+t_num), :, :) = temp_sozai;
end
 % 此处是否应该添加一个Hilbert包络线操作

% 此时sozai构建已经完成，开始对dim2进行permutation test

% 3种模式
% mod1:Rest vs. Imag    12    1
% mod2:Imag vs. Exec    23    2
% mod3:Rest vs. Exec    13    3

kekka = zeros(length(use_channels),length(Bandrange),3);    % 输出p值检验结果

for chn = 1:size(sozai,1)
    for band = 1:length(Bandrange)
        for m = 1:3
            id1 = m;
            id2 = mod(m,3)+1;
            seq1 = sozai(chn,:,band,id1);
            seq2 = sozai(chn,:,band,id2);         % 2段测试样本
            l1 = length(seq1);
            l2 = length(seq2); 
            label = [zeros(1,l1),ones(1,l2)];
            r_obs = corr([seq1 seq2]', label', 'type','spearman');
            tmp = [seq1 seq2];
            r_pdf = zeros(1,number_of_repetitions);     
            for i = 1:number_of_repetitions
                w = tmp(randperm(length(tmp))); 
                r_pdf(i) = corr(w', label', 'type','spearman'); 
            end
            p_value = 2*normcdf(abs(r_obs), mean(r_pdf, 1), std(r_pdf, 0, 1));
            kekka(chn,band,m) = p_value;

        end
    end
end



