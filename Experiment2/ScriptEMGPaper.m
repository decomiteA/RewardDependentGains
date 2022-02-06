% This script plots figure 5 and computes all the associated statistics of
% the following paper :
% Reward-dependent selection of feedback gains impacts rapid motor decision
% De Comite A., Crevecoeur F., Lefèvre P.
%
% Author : Antoine DE COMITE
% antoinedecomite@gmail.com
%
%
close all; clc; clear all;

%% 0. Loading and preprocessing the data
clc; clear all; close all;

DataM2Bis = load('PilotesBis.mat');
nSubjects = length(fieldnames(DataM2Bis));
Subjects = cell(nSubjects,1);
TotSubjects = zeros(480*nSubjects,1);
for ii = 1 : nSubjects
    Subjects{ii} = strcat('S',num2str(ii));
    TotSubjects((ii-1)*480+1:ii*480) = ii*ones(480,1);
end

mappingTP = [1 2 4 5 7 8 10 11 13 14];
TPright = [4 5 13 14];
TPleft = [7 8 10 11];
TotVectorTP = zeros(nSubjects*480,1);
TotReached = zeros(nSubjects*480,1);
TotSuccess = zeros(nSubjects*480,1);
TotReward = zeros(nSubjects*480,1);
TotForce = zeros(nSubjects*480,1);
TotRight = zeros(nSubjects*480,1);
for ii = 1 : nSubjects
    TotReached((ii-1)*480+1:ii*480) = BooleanReachedM2Bis(DataM2Bis.(Subjects{ii}));
    TotSuccess((ii-1)*480+1:ii*480) = BooleanSuccessM2Bis(DataM2Bis.(Subjects{ii}));
    TotVectorTP((ii-1)*480+1:ii*480) = DataM2Bis.(Subjects{ii}).vector_TP;
    TotIdem((ii-1)*480+1:ii*480) = mod(TotVectorTP((ii-1)*480+1:ii*480),3)-1;
    TotRight((ii-1)*480+1:ii*480) = 1*(TotVectorTP((ii-1)*480+1:ii*480)>12 | TotVectorTP((ii-1)*480+1:ii*480)<6) + 2*(TotVectorTP((ii-1)*480+1:ii*480)<12 & TotVectorTP((ii-1)*480+1:ii*480)>6);
    TotRight((ii-1)*480+1:ii*480) = (TotVectorTP((ii-1)*480+1:ii*480)>3).*TotRight((ii-1)*480+1:ii*480);
end

for ii = 1 : length(TotForce)
    switch TotVectorTP(ii)
        case {1,2}
            TotForce(ii) = 0;
        case {4,5}
            TotForce(ii) = 2;
        case {7,8}
            TotForce(ii) = -2;
        case {10,11}
            TotForce(ii) = 1;
        case {13,14}
            TotForce(ii) = -1;
    end
end

extremeMin = DataM2Bis.(Subjects{1}).time_matrix(1);
extremeMax = DataM2Bis.(Subjects{1}).time_matrix(end);
for ii = 2 : nSubjects
    extremeMin = max(extremeMin, DataM2Bis.(Subjects{ii}).time_matrix(1));
    extremeMax = min(extremeMax, DataM2Bis.(Subjects{ii}).time_matrix(end));
end
TimeVector = extremeMin : extremeMax;
TotEMG = zeros(480*nSubjects, length(TimeVector),2);
TotTimeMat = zeros(480*nSubjects,size(DataM2Bis.S1.matrix_timing,2));
for ii = 1 : nSubjects
    TotTimeMat((ii-1)*480+1:ii*480,:) = DataM2Bis.(Subjects{ii}).matrix_timing;
    TotEMG((ii-1)*480+1:ii*480,:,:) = DataM2Bis.(Subjects{ii}).matrix_EMG(:,find(DataM2Bis.(Subjects{ii}).time_matrix==TimeVector(1)):find(DataM2Bis.(Subjects{ii}).time_matrix==TimeVector(end)),:);
end
TotNearestTarget = zeros(480*nSubjects,1);
pos_targets = [-0.01 0.28; 0.09 0.28; 0.19 0.28];
for ii = 1 : nSubjects
    vectrick = DataM2Bis.(Subjects{ii}).time_matrix;
    for jj = 1 : 480
        if (TotReached((ii-1)*480+jj)==1)
            idxend = find(DataM2Bis.(Subjects{ii}).time_matrix == (min(abs(floor((DataM2Bis.(Subjects{ii}).matrix_timing(jj,8))*1000 - (DataM2Bis.(Subjects{ii}).matrix_timing(jj,6))*1000)),floor(vectrick(end)))));
            TotNearestTarget((ii-1)*480+jj) = NearestTarget([DataM2Bis.(Subjects{ii}).matrix_kinematics(jj,idxend,1), DataM2Bis.(Subjects{ii}).matrix_kinematics(jj,idxend,2)],pos_targets);
        end
    end
end

TotSwitch = ~(TotNearestTarget == 2 | TotNearestTarget ==0);

TotKine = zeros(480*nSubjects,length(TimeVector),2);
for ii = 1 : nSubjects
    TotKine((ii-1)*480+1:ii*480,:,:) = DataM2Bis.(Subjects{ii}).matrix_kinematics(:,find(DataM2Bis.(Subjects{ii}).time_matrix==TimeVector(1)):find(DataM2Bis.(Subjects{ii}).time_matrix==TimeVector(end)),:);
end

TotSpeed = zeros(480*nSubjects,length(TimeVector),2);
TotSpeed(:,3:end-2,:) = (TotKine(:,1:end-4,:)-8*TotKine(:,2:end-3,:)+8*TotKine(:,4:end-1,:)-TotKine(:,5:end,:))/(12*0.001);
TotReward = mod(TotVectorTP,3)-1;


%% Stats for pre activation of muscles

table48 = table(mean(TotEMG((TotVectorTP==mappingTP(3) | TotVectorTP==mappingTP(7)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),TotSubjects((TotVectorTP==mappingTP(3) | TotVectorTP==mappingTP(7)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)), TotNearestTarget((TotVectorTP==mappingTP(3) | TotVectorTP==mappingTP(7)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme48 = fitlme(table48,'EMG~Target+(1|Subjects)');

table610 = table(mean(TotEMG((TotVectorTP==mappingTP(5) | TotVectorTP==mappingTP(9)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),TotSubjects((TotVectorTP==mappingTP(5) | TotVectorTP==mappingTP(9)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)), TotNearestTarget((TotVectorTP==mappingTP(5) | TotVectorTP==mappingTP(9)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme610 = fitlme(table610,'EMG~Target+(1|Subjects)');

% The same with force levels 
table_large_del = table(mean(TotEMG((TotVectorTP==mappingTP(5)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),TotSubjects((TotVectorTP==mappingTP(5)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget((TotVectorTP==mappingTP(5)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Sujets','Condition'});
lme_large_del = fitlme(table_large_del,'EMG~Condition+(1|Sujets)');


table59 = table(mean(TotEMG((TotVectorTP==mappingTP(4) | TotVectorTP==mappingTP(8)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),TotSubjects((TotVectorTP==mappingTP(4) | TotVectorTP==mappingTP(8)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)), TotNearestTarget((TotVectorTP==mappingTP(4) | TotVectorTP==mappingTP(8)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme59 = fitlme(table59,'EMG~Target+(1|Subjects)');

table711 = table(mean(TotEMG((TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),TotSubjects((TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)), TotNearestTarget((TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme711 = fitlme(table711,'EMG~Target+(1|Subjects)');


%% Mean and SEM emg at movement onset 
figure('Name','Mean EMG at movement onset ','units','normalized'); hold on;
subplot(6,12,[52 53 54 64 65 66]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');
plot([1,2], [mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),1),mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==0),2),1)],'k.-','LineWidth',3,'MarkerSize',35);
plot([1,1], [mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),1)+std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),0,1)/sqrt(nSubjects),mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),1)-std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),0,1)/sqrt(nSubjects)],'k','LineWidth',3);
plot([2,2], [mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==0),2),1)+std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==0),2),0,1)/sqrt(nSubjects),mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==0),2),1)-std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),0,1)/sqrt(nSubjects)],'k','LineWidth',3);
xlim([0.7 2.3]); ylim([0.35 0.45]);

subplot(6,12,[55 56 57 67 68 69]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YColor','none');
plot([1,2], [mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),1),mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),1)],'k.-','LineWidth',3,'MarkerSize',35);
plot([1,1], [mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),1)+std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),0,1)/sqrt(nSubjects),mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),1)-std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),0,1)/sqrt(nSubjects)],'k','LineWidth',3);
plot([2,2], [mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),1)+std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),0,1)/sqrt(nSubjects),mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),1)-std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),0,1)/sqrt(nSubjects)],'k','LineWidth',3);
xlim([0.7 2.3]); ylim([0.35 0.45]);


%% Mean and SEM velocity at the beginning of movement
figure('Name','Velocity at perturbation onset'); hold on;
subplot(1,2,1); hold on; set(gca,'LineWidth',3); set(gca,'FontSize',14); set(gca,'Color','none'); hold on; xlim([-200 50]); ylim([-0.05 0.30]);
plot(-200:300,mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==300),2),1) - mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==-200):find(TimeVector==300),2),1),'k','Linewidth',3);
meanstd1 = mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==-200):find(TimeVector==300),2),1);
plot(-200:300,mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==300),2),1)-meanstd1 + std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==300),2)-meanstd1)/sqrt(nSubjects),'k:','LineWidth',3)
plot(-200:300,mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==300),2),1)-meanstd1 - std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==300),2)-meanstd1)/sqrt(nSubjects),'k:','LineWidth',3)
yline(0,'k','LineWidth',3); ylabel('Forward velocity [cm/s]');


subplot(1,2,2); hold on; set(gca,'LineWidth',3); set(gca,'FontSize',14); set(gca,'Color','none'); hold on; xlim([-200 50]); ylim([-0.05 0.30]);
plot(-200:300,mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==300),2),1) - mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==-200):find(TimeVector==300),2),1),'k','Linewidth',3);
meanstd2 = mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==-200):find(TimeVector==300),2),1);
plot(-200:300,mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==300),2),1)-meanstd2 + std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==300),2)-meanstd2)/sqrt(nSubjects),'k:','LineWidth',3)
plot(-200:300,mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==300),2),1)-meanstd2 - std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==300),2)-meanstd2)/sqrt(nSubjects),'k:','LineWidth',3)
yline(0,'k','LineWidth',3);

%%

tablespeed48 = table(TotSpeed((TotVectorTP==mappingTP(4)|TotVectorTP==mappingTP(8)) & (TotNearestTarget==2 | TotNearestTarget==3),TimeVector==0,2),TotNearestTarget((TotVectorTP==mappingTP(4)|TotVectorTP==mappingTP(8)) & (TotNearestTarget==2 | TotNearestTarget==3)),TotSubjects((TotVectorTP==mappingTP(4)|TotVectorTP==mappingTP(8)) & (TotNearestTarget==2 | TotNearestTarget==3)),'VariableNames',{'Speed','TP','Subjects'});
lmespeed48 = fitlme(tablespeed48,'Speed~TP+(1|Subjects)');

tablespeed610 = table(TotSpeed((TotVectorTP==mappingTP(6)|TotVectorTP==mappingTP(10)) & (TotNearestTarget==2 | TotNearestTarget==1),TimeVector==0,2),TotNearestTarget((TotVectorTP==mappingTP(6)|TotVectorTP==mappingTP(10)) & (TotNearestTarget==2 | TotNearestTarget==1)),TotSubjects((TotVectorTP==mappingTP(6)|TotVectorTP==mappingTP(10)) & (TotNearestTarget==2 | TotNearestTarget==1)),'VariableNames',{'Speed','TP','Subjects'});
lmespeed610 = fitlme(tablespeed610,'Speed~TP+(1|Subjects)');

tablespeed59 = table(TotSpeed((TotVectorTP==mappingTP(3)|TotVectorTP==mappingTP(7)) & (TotNearestTarget==2 | TotNearestTarget==3),TimeVector==0,2),TotNearestTarget((TotVectorTP==mappingTP(3)|TotVectorTP==mappingTP(7)) & (TotNearestTarget==2 | TotNearestTarget==3)),TotSubjects((TotVectorTP==mappingTP(3)|TotVectorTP==mappingTP(7)) & (TotNearestTarget==2 | TotNearestTarget==3)),'VariableNames',{'Speed','TP','Subjects'});
lmespeed59 = fitlme(tablespeed59,'Speed~TP+(1|Subjects)');

tablespeed711 = table(TotSpeed((TotVectorTP==mappingTP(5)|TotVectorTP==mappingTP(9)) & (TotNearestTarget==2 | TotNearestTarget==1),TimeVector==0,2),TotNearestTarget((TotVectorTP==mappingTP(5)|TotVectorTP==mappingTP(9)) & (TotNearestTarget==2 | TotNearestTarget==1)),TotSubjects((TotVectorTP==mappingTP(5)|TotVectorTP==mappingTP(9)) & (TotNearestTarget==2 | TotNearestTarget==1)),'VariableNames',{'Speed','TP','Subjects'});
lmespeed711 = fitlme(tablespeed711,'Speed~TP+(1|Subjects)');


tablekine48 = table(TotKine((TotVectorTP==mappingTP(4)|TotVectorTP==mappingTP(8)) & (TotNearestTarget==2 | TotNearestTarget==3),TimeVector==0,1),TotNearestTarget((TotVectorTP==mappingTP(4)|TotVectorTP==mappingTP(8)) & (TotNearestTarget==2 | TotNearestTarget==3)),TotSubjects((TotVectorTP==mappingTP(4)|TotVectorTP==mappingTP(8)) & (TotNearestTarget==2 | TotNearestTarget==3)),'VariableNames',{'Pos','TP','Subjects'});
lmekine48 = fitlme(tablekine48,'Pos~TP+(1|Subjects)');

tablekine610 = table(TotKine((TotVectorTP==mappingTP(6)|TotVectorTP==mappingTP(10)) & (TotNearestTarget==2 | TotNearestTarget==1),TimeVector==0,1),TotNearestTarget((TotVectorTP==mappingTP(6)|TotVectorTP==mappingTP(10)) & (TotNearestTarget==2 | TotNearestTarget==1)),TotSubjects((TotVectorTP==mappingTP(6)|TotVectorTP==mappingTP(10)) & (TotNearestTarget==2 | TotNearestTarget==1)),'VariableNames',{'Pos','TP','Subjects'});
lmekine610 = fitlme(tablekine610,'Pos~TP+(1|Subjects)');
%% EMG activity binned before movement initiation
figure('Name','Binned EMG traces');
subplot(6,12,[9 10 21 22]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');
bar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1)],'ShowBaseLine','off','FaceColor',[1 0 0],'EdgeColor','none');  
xlim([0.3 2.7]); set(gca,'YColor','none');
er = errorbar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1)],[std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2)),std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2))]/sqrt(nSubjects),[std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2)),std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2))]/sqrt(nSubjects));
er.Color = [0 0 0]; er.LineWidth=3;

table48pec = table(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & (TotNearestTarget==2 | TotNearestTarget==3) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),TotSubjects((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & (TotNearestTarget==2 | TotNearestTarget==3) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & (TotNearestTarget==2 | TotNearestTarget==3) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme48pec = fitlme(table48pec,'EMG~Target+(1|Subjects)');

subplot(6,12,[11 12 23 24]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');
bar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1)],'ShowBaseLine','off','FaceColor',[1 0 0],'EdgeColor','none');  
xlim([0.3 2.7]); set(gca,'YAxisLocation','right');ylabel('EMG activity [a.u.]');
er = errorbar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1)],[std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2)),std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2))]/sqrt(nSubjects),[std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2)),std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2))]/sqrt(nSubjects));
er.Color = [0 0 0]; er.LineWidth=3;

table610pec = table(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),TotSubjects((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme610pec = fitlme(table610pec,'EMG~Target+(1|Subjects)');


subplot(6,12,[33 34 45 46]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');
bar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1)],'ShowBaseLine','off','FaceColor',[1 0 0],'EdgeColor','none');  
xlim([0.3 2.7]); set(gca,'YColor','none'); xticks([1 2]); xticklabels({'Center','Lateral'}); xtickangle(90);
er = errorbar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1)],[std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2)),std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2))]/sqrt(nSubjects),[std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2)),std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2))]/sqrt(nSubjects));
er.Color = [0 0 0]; er.LineWidth=3;

table48del = table(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & (TotNearestTarget==2 | TotNearestTarget==3) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),TotSubjects((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & (TotNearestTarget==2 | TotNearestTarget==3) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & (TotNearestTarget==2 | TotNearestTarget==3) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme48del = fitlme(table48del,'EMG~Target+(1|Subjects)');

subplot(6,12,[35 36 47 48]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');
bar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1)],'ShowBaseLine','off','FaceColor',[1 0 0],'EdgeColor','none');  
xlim([0.3 2.7]); set(gca,'YAxisLocation','right'); xticks([1 2]); xticklabels({'Center','Lateral'}); xtickangle(90);
er = errorbar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1)],[std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2)),std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2))]/sqrt(nSubjects),[std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2)),std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2))]/sqrt(nSubjects));
er.Color = [0 0 0]; er.LineWidth=3;

table610del = table(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),TotSubjects((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme610del = fitlme(table610del,'EMG~Target+(1|Subjects)');

figure('Name','Representation of mean EMG traces','units','normalized','outerposition',[0 0 1 1]); hold on;
subplot(2,2,1); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:500, movmean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==500),1),1),1),'b','LineWidth',4);
plot(-200:500, movmean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==500),1),1),1),'g','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 1.5]); yLim = ylim;

subplot(2,2,2); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:500, movmean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==500),1),1),1),'b','LineWidth',4);
plot(-200:500, movmean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==500),1),1),1),'r','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 1.5]);yLim = ylim;

subplot(2,2,3); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:500, movmean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==500),2),1),1),'b','LineWidth',4);
plot(-200:500, movmean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==500),2),1),1),'g','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 1.5]);yLim = ylim;

subplot(2,2,4); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:500, movmean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==500),2),1),1),'b','LineWidth',4);
plot(-200:500, movmean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==500),2),1),1),'r','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 1.5]);yLim = ylim;