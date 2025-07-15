function MakeEMGplot(pn,carind)
filepath=['E:\Schalk\Data_Analysis\P',num2str(pn),'\Data_Process_Standard\2_Data_Alignment_Resize\'];
datasetname=strcat('Power_P',num2str(pn),'_H1_S12_C12345_B60_140.mat');
if carind==0
    load([filepath,'WITHOUTCAR\',datasetname]);
elseif carind==1
    load([filepath,'WITHCAR\',datasetname]);
elseif carind==2
    load([filepath,'WITHCAR_CW_GROUPED\',datasetname]);
elseif carind==3
    load([filepath,'LOCALREF\',datasetname]);
elseif carind==4
    load([filepath,'WITHCAR_Median\',datasetname]);
elseif carind==5
    load([filepath,'WITHCAR_Chn_Spec\',datasetname]);
elseif carind==6
    load([filepath,'WITHCAR_Ele_Spec\',datasetname]);
elseif carind==7
    load([filepath,'BiPolar\',datasetname]);
else
end
EMGOnset=Power_Data(:,:,size(Power_Data,3)-1);
fealabel=Power_Data(:,:,size(Power_Data,3));
EMGinex=[];
emg_trial=Power_Data(:,:,size(Power_Data,3)-2);
for i=1:size(EMGOnset,2)
    if isempty(find(EMGOnset(:,i)==1))
        EMGinex(i)=(560+2000)/5;
    else
        EMGinex(i)=find(EMGOnset(:,i)==1);
    end
fealabelind(i)=min(find(fealabel(:,i)>0));
end
delay=EMGinex-fealabelind;
interv = 1;
trialNo = size(emg_trial,2);
IC = zeros(trialNo,1);
for i=1:max(max(fealabel))
    index = (i-1)*20+1:i*20;
    [delay(index),IC(index)] = sort(delay(index),'ascend');
    IC(index) = IC(index) + (i-1)*20;
end
fealabelind=fealabelind(IC);
emg_trial = emg_trial(:,IC);
emg_trial = (mapminmax(emg_trial',0,1))';
delay(delay<=0) = 1;
filename=['S',num2str(pn),'_EMGonsetVerification.jpg'];
figurepath=[]
HF=figure; %set(gca,'position',[0,0 1,1])
axes('position',[0.05 0.05 0.94 0.94])
plot(emg_trial+repmat((0:trialNo-1)*interv,size(emg_trial,1),1),'k');axis tight
hold on;
dy = arrayfun(@(i) emg_trial(delay(i),i),1:trialNo);
for j=1:trialNo
plot([delay(j)+fealabelind(j) delay(j)+fealabelind(j)],[dy(j)+(j-1)*interv-0.5 dy(j)+(j-1)*interv+0.5],'r','linewidth',1.5);
end
print(HF,'-djpeg','-r300',['E:\Schalk\Data_Analysis\Power_Transfer_Ave\GOSDETCT\Figures\EMG_ONSET\',filename]);
close(HF);
end