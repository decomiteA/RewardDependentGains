% Script for the EMG analysis in the paper

% Antoine De Comite
% 24th of March 2020


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
%% Let's take a look at these EMG
TotEMGIndMean = ComputeIndMeanCdt(TotEMG,TotVectorTP,TotNearestTarget,TotSubjects,nSubjects);
SubjectToPlot=[1 2 3 4 6 7 8 9 11 12 13 14 15 17 19]; % the other data was corrupted or not exploitable


% Figure for the paper
close all;
figure('Name','Graphical representation of EMG correlates','units','normalized','outerposition',[0 0 1 1]); hold on;
subplot(2,2,1); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none'); title('Small mechanical perturbation');
plot(-200:500, movmean(customMean(TotEMGIndMean(8,find(TimeVector==-200):find(TimeVector==500),2,1,SubjectToPlot),5,0),5),'b','LineWidth',4);
plot(-200:500, movmean(customMean(TotEMGIndMean(8,find(TimeVector==-200):find(TimeVector==500),3,1,SubjectToPlot),5,0),5),'g','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 1.5]); yLim = ylim;
plot(linspace(0,0,100),linspace(yLim(1),yLim(end),100),'k:','LineWidth',2.5);
ylabel('EMG activity [a.u.]');


subplot(2,2,2); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none'); title('Large mechanical perturbation');
plot(-200:500, movmean(customMean(TotEMGIndMean(4,find(TimeVector==-200):find(TimeVector==500),2,1,SubjectToPlot),5,0),5),'b','LineWidth',4);
plot(-200:500, movmean(customMean(TotEMGIndMean(4,find(TimeVector==-200):find(TimeVector==500),3,1,SubjectToPlot),5,0),5),'g','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 1.5]);yLim = ylim;
ylabel('EMG activity [a.u.]');plot(linspace(0,0,100),linspace(yLim(1),yLim(end),100),'k:','LineWidth',2.5);

subplot(2,2,3); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:500, movmean(customMean(TotEMGIndMean(10,find(TimeVector==-200):find(TimeVector==500),2,2,SubjectToPlot),5,0),5),'b','LineWidth',4);
plot(-200:500, movmean(customMean(TotEMGIndMean(10,find(TimeVector==-200):find(TimeVector==500),1,2,SubjectToPlot),5,0),5),'r','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 2]);yLim = ylim;
ylabel('EMG activity [a.u.]');plot(linspace(0,0,100),linspace(yLim(1),yLim(end),100),'k:','LineWidth',2.5);

subplot(2,2,4); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:500, movmean(customMean(TotEMGIndMean(6,find(TimeVector==-200):find(TimeVector==500),2,2,SubjectToPlot),5,0),5),'b','LineWidth',4);
plot(-200:500, movmean(customMean(TotEMGIndMean(6,find(TimeVector==-200):find(TimeVector==500),1,2,SubjectToPlot),5,0),5),'r','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 2]);yLim = ylim;
ylabel('EMG activity [a.u.]');plot(linspace(0,0,100),linspace(yLim(1),yLim(end),100),'k:','LineWidth',2.5);




figure('Name','Modification post round 1','units','normalized','outerposition',[0 0 1 1]); hold on;
subplot(2,2,1); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:500, movmean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-200):find(TimeVector==500),2,1,SubjectToPlot),5,0),5,0),8),'b','LineWidth',4);
plot(-200:500, movmean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-200):find(TimeVector==500),3,1,SubjectToPlot),5,0),5,0),8),'g','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 1.5]); yLim = ylim;

subplot(2,2,2); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:500, movmean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-200):find(TimeVector==500),2,1,SubjectToPlot),5,0),5,0),8),'b','LineWidth',4);
plot(-200:500, movmean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-200):find(TimeVector==500),1,1,SubjectToPlot),5,0),5,0),8),'r','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 1.5]);yLim = ylim;

subplot(2,2,3); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:500, movmean(customMean(customMean(TotEMGIndMean([4 8],find(TimeVector==-200):find(TimeVector==500),2,2,SubjectToPlot),5,0),5,0),8),'b','LineWidth',4);
plot(-200:500, movmean(customMean(customMean(TotEMGIndMean([4 8],find(TimeVector==-200):find(TimeVector==500),3,2,SubjectToPlot),5,0),5,0),8),'g','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 1.5]);yLim = ylim;

subplot(2,2,4); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:500, movmean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-200):find(TimeVector==500),2,2,SubjectToPlot),5,0),5,0),8),'b','LineWidth',4);
plot(-200:500, movmean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-200):find(TimeVector==500),1,2,SubjectToPlot),5,0),5,0),8),'r','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-150','0','150','300','450'}); xlabel('Time [ms]'); ylim([0 1.5]);yLim = ylim;


figure('Name','Test review question','units','normalized');
subplot(2,2,1); hold on; set(gca,'Color','none'); set(gca,'LineWidth',2); set(gca,'FontSize',14);
plot([1,2],[mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),2,1,SubjectToPlot),5,0),5,0),2),1),mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),3,1,SubjectToPlot),5,0),5,0),2),1)],'k.-','MarkerSize',35,'LineWidth',2);
plot([1,1],[mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),2,1,SubjectToPlot),5,0),5,0),2),1)+mean(std(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),2,1,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects),mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),2,1,SubjectToPlot),5,0),5,0),2),1)-mean(std(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),2,1,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects)],'k-','LineWidth',2);
plot([2,2],[mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),3,1,SubjectToPlot),5,0),5,0),2),1)+mean(std(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),3,1,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects),mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),3,1,SubjectToPlot),5,0),5,0),2),1)-mean(std(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),3,1,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects)],'k-','LineWidth',2);
xlim([0 3]);


subplot(2,2,3); hold on; set(gca,'Color','none'); set(gca,'LineWidth',2); set(gca,'FontSize',14);
plot([1,2],[mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),2,1,SubjectToPlot),5,0),5,0),2),1),mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),1,1,SubjectToPlot),5,0),5,0),2),1)],'k.-','MarkerSize',35,'LineWidth',2);
plot([1,1],[mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),2,1,SubjectToPlot),5,0),5,0),2),1)+mean(std(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),2,1,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects),mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),2,1,SubjectToPlot),5,0),5,0),2),1)-mean(std(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),2,1,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects)],'k-','LineWidth',2);
plot([2,2],[mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),1,1,SubjectToPlot),5,0),5,0),2),1)+mean(std(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),1,1,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects),mean(mean(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),1,1,SubjectToPlot),5,0),5,0),2),1)-mean(std(customMean(customMean(TotEMGIndMean([8 4],find(TimeVector==-150):find(TimeVector==0),1,1,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects)],'k-','LineWidth',2);
xlim([0 3]);

subplot(2,2,2); hold on; set(gca,'Color','none'); set(gca,'LineWidth',2); set(gca,'FontSize',14);
plot([1,2],[mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),2,2,SubjectToPlot),5,0),5,0),2),1),mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),3,2,SubjectToPlot),5,0),5,0),2),1)],'k.-','MarkerSize',35,'LineWidth',2);
plot([1,1],[mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),2,2,SubjectToPlot),5,0),5,0),2),1)+mean(std(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),2,2,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects),mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),2,2,SubjectToPlot),5,0),5,0),2),1)-mean(std(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),2,2,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects)],'k-','LineWidth',2);
plot([2,2],[mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),3,2,SubjectToPlot),5,0),5,0),2),1)+mean(std(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),3,2,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects),mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),3,2,SubjectToPlot),5,0),5,0),2),1)-mean(std(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),3,2,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects)],'k-','LineWidth',2);
xlim([0 3]);

subplot(2,2,4); hold on; set(gca,'Color','none'); set(gca,'LineWidth',2); set(gca,'FontSize',14);
plot([1,2],[mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),2,2,SubjectToPlot),5,0),5,0),2),1),mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),1,2,SubjectToPlot),5,0),5,0),2),1)],'k.-','MarkerSize',35,'LineWidth',2);
plot([1,1],[mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),2,2,SubjectToPlot),5,0),5,0),2),1)+mean(std(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),2,2,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects),mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),2,2,SubjectToPlot),5,0),5,0),2),1)-mean(std(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),2,2,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects)],'k-','LineWidth',2);
plot([2,2],[mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),1,2,SubjectToPlot),5,0),5,0),2),1)+mean(std(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),1,2,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects),mean(mean(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),1,2,SubjectToPlot),5,0),5,0),2),1)-mean(std(customMean(customMean(TotEMGIndMean([6 10],find(TimeVector==-150):find(TimeVector==0),1,2,SubjectToPlot),5,0),5,0),0,2),1)/sqrt(nSubjects)],'k-','LineWidth',2);
xlim([0 3]);
%% Investigation of what happen before the force in the EMG data, let's look at this for all subjects (except for those that didn't show switch) - has to be robust


% Cocontraction - Time bins

TimeBins = [-150 -100;
    -100 -50 ;
    -50    0];

BinnedEMGValues = zeros(size(TotEMGIndMean,1),size(TimeBins,1),3,2,nSubjects);% 1. Cdt, 2. time bin, 3. target, 4. muscles, 5. Sujets
for ii = 1 : size(TimeBins,1)
    BinnedEMGValues(:,ii,:,:,:) = mean(TotEMGIndMean(:,find(TimeVector==TimeBins(ii,1)):find(TimeVector==TimeBins(ii,2)),:,:,:),2);
end


figure('Name','Investigation the binned EMG value','units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot([1 2 3], customMean(BinnedEMGValues(8,:,2,1,SubjectToPlot),5,0),'.b','MarkerSize',25);
plot([1.05 2.05 3.05], customMean(BinnedEMGValues(8,:,3,1,SubjectToPlot),5,0),'.g','MarkerSize',25);
for ii = 1 : size(BinnedEMGValues,2)
    plot(linspace(ii,ii), linspace(customMean(BinnedEMGValues(8,ii,2,1,SubjectToPlot),5,0)+customSTD(BinnedEMGValues(8,ii,2,1,SubjectToPlot),5,0)/sqrt(nSubjects),customMean(BinnedEMGValues(8,ii,2,1,SubjectToPlot),5,0)-customSTD(BinnedEMGValues(8,ii,2,1,SubjectToPlot),5,0)/sqrt(nSubjects)),'b','LineWidth',2);
    plot(linspace(ii+0.05,ii+0.05), linspace(customMean(BinnedEMGValues(8,ii,3,1,SubjectToPlot),5,0)+customSTD(BinnedEMGValues(8,ii,3,1,SubjectToPlot),5,0)/sqrt(nSubjects),customMean(BinnedEMGValues(8,ii,3,1,SubjectToPlot),5,0)-customSTD(BinnedEMGValues(8,ii,3,1,SubjectToPlot),5,0)/sqrt(nSubjects)),'g','LineWidth',2);
end
xlim([0.7 3.3]); ylim([0 0.5]);

subplot(2,2,2); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot([1 2 3], customMean(BinnedEMGValues(4,:,2,1,SubjectToPlot),5,0),'.b','MarkerSize',25);
plot([1.05 2.05 3.05], customMean(BinnedEMGValues(4,:,3,1,SubjectToPlot),5,0),'.g','MarkerSize',25);
for ii = 1 : size(BinnedEMGValues,2)
    plot(linspace(ii,ii), linspace(customMean(BinnedEMGValues(4,ii,2,1,SubjectToPlot),5,0)+customSTD(BinnedEMGValues(4,ii,2,1,SubjectToPlot),5,0)/sqrt(nSubjects),customMean(BinnedEMGValues(4,ii,2,1,SubjectToPlot),5,0)-customSTD(BinnedEMGValues(4,ii,2,1,SubjectToPlot),5,0)/sqrt(nSubjects)),'b','LineWidth',2);
    plot(linspace(ii+0.05,ii+0.05), linspace(customMean(BinnedEMGValues(4,ii,3,1,SubjectToPlot),5,0)+customSTD(BinnedEMGValues(4,ii,3,1,SubjectToPlot),5,0)/sqrt(nSubjects),customMean(BinnedEMGValues(4,ii,3,1,SubjectToPlot),5,0)-customSTD(BinnedEMGValues(4,ii,3,1,SubjectToPlot),5,0)/sqrt(nSubjects)),'g','LineWidth',2);
end
xlim([0.7 3.3]); ylim([0 0.5]);

subplot(2,2,3); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot([1 2 3], customMean(BinnedEMGValues(10,:,2,2,SubjectToPlot),5,0),'.b','MarkerSize',25);
plot([1.05 2.05 3.05], customMean(BinnedEMGValues(10,:,1,2,SubjectToPlot),5,0),'.r','MarkerSize',25);
for ii = 1 : size(BinnedEMGValues,2)
    plot(linspace(ii,ii), linspace(customMean(BinnedEMGValues(10,ii,2,2,SubjectToPlot),5,0)+customSTD(BinnedEMGValues(10,ii,2,2,SubjectToPlot),5,0)/sqrt(nSubjects),customMean(BinnedEMGValues(10,ii,2,2,SubjectToPlot),5,0)-customSTD(BinnedEMGValues(10,ii,2,2,SubjectToPlot),5,0)/sqrt(nSubjects)),'b','LineWidth',2);
    plot(linspace(ii+0.05,ii+0.05), linspace(customMean(BinnedEMGValues(10,ii,1,2,SubjectToPlot),5,0)+customSTD(BinnedEMGValues(10,ii,1,2,SubjectToPlot),5,0)/sqrt(nSubjects),customMean(BinnedEMGValues(10,ii,1,2,SubjectToPlot),5,0)-customSTD(BinnedEMGValues(10,ii,1,2,SubjectToPlot),5,0)/sqrt(nSubjects)),'r','LineWidth',2);
end
xlim([0.7 3.3]); ylim([0 0.5]);

subplot(2,2,4); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot([1 2 3], customMean(BinnedEMGValues(6,:,2,2,SubjectToPlot),5,0),'.b','MarkerSize',25);
plot([1.05 2.05 3.05], customMean(BinnedEMGValues(6,:,1,2,SubjectToPlot),5,0),'.r','MarkerSize',25);
for ii = 1 : size(BinnedEMGValues,2)
    plot(linspace(ii,ii), linspace(customMean(BinnedEMGValues(6,ii,2,2,SubjectToPlot),5,0)+customSTD(BinnedEMGValues(6,ii,2,2,SubjectToPlot),5,0)/sqrt(nSubjects),customMean(BinnedEMGValues(6,ii,2,2,SubjectToPlot),5,0)-customSTD(BinnedEMGValues(6,ii,2,2,SubjectToPlot),5,0)/sqrt(nSubjects)),'b','LineWidth',2);
    plot(linspace(ii+0.05,ii+0.05), linspace(customMean(BinnedEMGValues(6,ii,1,2,SubjectToPlot),5,0)+customSTD(BinnedEMGValues(6,ii,1,2,SubjectToPlot),5,0)/sqrt(nSubjects),customMean(BinnedEMGValues(6,ii,1,2,SubjectToPlot),5,0)-customSTD(BinnedEMGValues(6,ii,1,2,SubjectToPlot),5,0)/sqrt(nSubjects)),'r','LineWidth',2);
end
xlim([0.7 3.3]); ylim([0 0.5]);





% Let's perform some statistical analysis
%%

% For trial number 8
% 5 10 16 18 20

tablecocontraction81 = table(mean(TotEMG(TotVectorTP==mappingTP(8) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==-100),1),2),TotSubjects(TotVectorTP==mappingTP(8) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(8) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme81 = fitlme(tablecocontraction81,'EMG~Target+(1|Subjects)');
tablecocontraction82 = table(mean(TotEMG(TotVectorTP==mappingTP(8) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-100):find(TimeVector==-50),1),2),TotSubjects(TotVectorTP==mappingTP(8) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(8) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme82 = fitlme(tablecocontraction82,'EMG~Target+(1|Subjects)');
tablecocontraction83 = table(mean(TotEMG(TotVectorTP==mappingTP(8) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-50):find(TimeVector==0),1),2),TotSubjects(TotVectorTP==mappingTP(8) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(8) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme83 = fitlme(tablecocontraction83,'EMG~Target+(1|Subjects)');


tablecocontraction41 = table(mean(TotEMG(TotVectorTP==mappingTP(4) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==-100),1),2),TotSubjects(TotVectorTP==mappingTP(4) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(4) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme41 = fitlme(tablecocontraction41,'EMG~Target+(1|Subjects)');
tablecocontraction42 = table(mean(TotEMG(TotVectorTP==mappingTP(4) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-100):find(TimeVector==-50),1),2),TotSubjects(TotVectorTP==mappingTP(4) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(4) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme42 = fitlme(tablecocontraction42,'EMG~Target+(1|Subjects)');
tablecocontraction43 = table(mean(TotEMG(TotVectorTP==mappingTP(4) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-50):find(TimeVector==0),1),2),TotSubjects(TotVectorTP==mappingTP(4) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(4) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme43 = fitlme(tablecocontraction43,'EMG~Target+(1|Subjects)');

tablecocontraction101 = table(mean(TotEMG(TotVectorTP==mappingTP(10) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==-100),2),2),TotSubjects(TotVectorTP==mappingTP(10) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(10) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme101 = fitlme(tablecocontraction101,'EMG~Target+(1|Subjects)');
tablecocontraction102 = table(mean(TotEMG(TotVectorTP==mappingTP(10) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-100):find(TimeVector==-50),2),2),TotSubjects(TotVectorTP==mappingTP(10) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(10) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme102 = fitlme(tablecocontraction102,'EMG~Target+(1|Subjects)');
tablecocontraction103 = table(mean(TotEMG(TotVectorTP==mappingTP(10) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-50):find(TimeVector==0),2),2),TotSubjects(TotVectorTP==mappingTP(10) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(10) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme103 = fitlme(tablecocontraction103,'EMG~Target+(1|Subjects)');


tablecocontraction61 = table(mean(TotEMG(TotVectorTP==mappingTP(6) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==-100),2),2),TotSubjects(TotVectorTP==mappingTP(6) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(6) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme61 = fitlme(tablecocontraction61,'EMG~Target+(1|Subjects)');
tablecocontraction62 = table(mean(TotEMG(TotVectorTP==mappingTP(6) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-100):find(TimeVector==-50),2),2),TotSubjects(TotVectorTP==mappingTP(6) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(6) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme62 = fitlme(tablecocontraction62,'EMG~Target+(1|Subjects)');
tablecocontraction63 = table(mean(TotEMG(TotVectorTP==mappingTP(6) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-50):find(TimeVector==0),2),2),TotSubjects(TotVectorTP==mappingTP(6) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget(TotVectorTP==mappingTP(6) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme63 = fitlme(tablecocontraction63,'EMG~Target+(1|Subjects)');

% Removing the "bad" subjects is good for significance
%% Final analysis (only 1 time bin)
% Let's now mixe large and small mechanical perturbations (since we only
% look at the co-contraction, we can do this because there was no clue of
% the intensity of the force yet). 

table48 = table(mean(TotEMG((TotVectorTP==mappingTP(3) | TotVectorTP==mappingTP(7)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),TotSubjects((TotVectorTP==mappingTP(3) | TotVectorTP==mappingTP(7)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)), TotNearestTarget((TotVectorTP==mappingTP(3) | TotVectorTP==mappingTP(7)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme48 = fitlme(table48,'EMG~Target+(1|Subjects)');

table610 = table(mean(TotEMG((TotVectorTP==mappingTP(5) | TotVectorTP==mappingTP(9)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),TotSubjects((TotVectorTP==mappingTP(5) | TotVectorTP==mappingTP(9)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)), TotNearestTarget((TotVectorTP==mappingTP(5) | TotVectorTP==mappingTP(9)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme610 = fitlme(table610,'EMG~Target+(1|Subjects)');

% The same with force levels 
table_large_del = table(mean(TotEMG((TotVectorTP==mappingTP(5)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),TotSubjects((TotVectorTP==mappingTP(5)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),TotNearestTarget((TotVectorTP==mappingTP(5)) & (TotNearestTarget==2 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Sujets','Condition'});
lme_large_del = fitlme(table_large_del,'EMG~Condition+(1|Sujets)');

%% Let's do this with the other condition

table59 = table(mean(TotEMG((TotVectorTP==mappingTP(4) | TotVectorTP==mappingTP(8)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),TotSubjects((TotVectorTP==mappingTP(4) | TotVectorTP==mappingTP(8)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)), TotNearestTarget((TotVectorTP==mappingTP(4) | TotVectorTP==mappingTP(8)) & (TotNearestTarget==3 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme59 = fitlme(table59,'EMG~Target+(1|Subjects)');

table711 = table(mean(TotEMG((TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),TotSubjects((TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)), TotNearestTarget((TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10)) & (TotNearestTarget==1 | TotNearestTarget==2) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20)),'VariableNames',{'EMG','Subjects','Target'});
lme711 = fitlme(table711,'EMG~Target+(1|Subjects)');


%% Figure for showing what's happening before movement
close all; 
figure('Name','Whats happening at the time of decision','units','normalized','outerposition',[0 0 1 1]);
subplot(6,12,[1 2 3 4 13 14 15 16]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');set(gca,'XColor','none');
plot(-200:500, movmean(customMean(TotEMGIndMean(8,find(TimeVector==-200):find(TimeVector==500),2,1,SubjectToPlot),5,0),5),'b','LineWidth',4);
plot(-200:500, movmean(customMean(TotEMGIndMean(8,find(TimeVector==-200):find(TimeVector==500),3,1,SubjectToPlot),5,0),5),'Color',[0 0.4 0],'LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-0.15','0','0.15','0.3','0.45'}); xlabel('Time [s]'); ylim([0 1.75]); yLim = ylim;
plot(linspace(0,0,100),linspace(yLim(1),yLim(end),100),'k:','LineWidth',2.5);
ylabel('EMG activity [a.u.]');


subplot(6,12,[25 26 27 28 37 38 39 40]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');
plot(-200:500, movmean(customMean(TotEMGIndMean(10,find(TimeVector==-200):find(TimeVector==500),2,2,SubjectToPlot),5,0),5),'b:','LineWidth',4);
plot(-200:500, movmean(customMean(TotEMGIndMean(10,find(TimeVector==-200):find(TimeVector==500),1,2,SubjectToPlot),5,0),5),'r','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-0.15','0','0.15','0.3','0.45'}); xlabel('Time [s]'); ylim([0 1.75]); yLim = ylim;
plot(linspace(0,0,100),linspace(yLim(1),yLim(end),100),'k:','LineWidth',2.5);


subplot(6,12,[5 6 7 8 17 18 19 20]);hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none'); set(gca,'YColor','none');
plot(-200:500, movmean(customMean(TotEMGIndMean(4,find(TimeVector==-200):find(TimeVector==500),2,1,SubjectToPlot),5,0),5),'b','LineWidth',4);
plot(-200:500, movmean(customMean(TotEMGIndMean(4,find(TimeVector==-200):find(TimeVector==500),3,1,SubjectToPlot),5,0),5),'Color',[0 0.4 0],'LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-0.15','0','0.15','0.3','0.45'}); xlabel('Time [s]'); ylim([0 1.75]); yLim = ylim;
plot(linspace(0,0,100),linspace(yLim(1),yLim(end),100),'k:','LineWidth',2.5);
ylabel('EMG activity [a.u.]');

subplot(6,12,[29 30 31 32 41 42 43 44]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YColor','none');
plot(-200:500, movmean(customMean(TotEMGIndMean(6,find(TimeVector==-200):find(TimeVector==500),2,2,SubjectToPlot),5,0),5),'b','LineWidth',4);
plot(-200:500, movmean(customMean(TotEMGIndMean(6,find(TimeVector==-200):find(TimeVector==500),1,2,SubjectToPlot),5,0),5),'r','LineWidth',4);
xlim([-200 500]); xticks([-150 0 150 300 450]); xticklabels({'-0.15','0','0.15','0.3','0.45'}); xlabel('Time [s]'); ylim([0 1.75]); yLim = ylim;
plot(linspace(0,0,100),linspace(yLim(1),yLim(end),100),'k:','LineWidth',2.5);
ylabel('EMG activity [a.u.]');
%%


% bar plots for the emg activity at perturbation onset 
% ajout des barres d'erreur

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Neuromatch talk !!! %%
%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
figure('Name','Graphical representation of the EMG trace','units','normalized','outerposition',[0 0 1 1]); hold on;
set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');

subplot(1,2,1); set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none'); hold on;
plot(-200:0,mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==0),1),1),'b','LineWidth',2.5);
plot(-200:0,mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(4)) & (TotNearestTarget==3 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==0),1),1),'r','LineWidth',2.5);
legend('Central','Lateral','FontSize',16);

subplot(1,2,2); set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none'); hold on;
plot(-200:0,mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==0),2),1),'b','LineWidth',2.5);
plot(-200:0,mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(4)) & (TotNearestTarget==3 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==0),2),1),'r','LineWidth',2.5);
legend('Central','Lateral','FontSize',16);

figure; hold on; 
plot(-200:500,mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==500),2),1),'b','LineWidth',2);
plot(-200:500,mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(4)) & (TotNearestTarget==3 | TotNearestTarget==1),find(TimeVector==-200):find(TimeVector==500),2),1),'r','LineWidth',2);
%%

figure('Name','Figure neuromatch 3.0','units','normalized','outerposition',[0.1 0.1 0.6 0.6]);
subplot(1,2,1); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:50,mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(4) | TotVectorTP==mappingTP(2)) & TotNearestTarget==2,find(TimeVector==-200):find(TimeVector==50),2),1),'k','LineWidth',2.5);
plot(-200:50,mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(4) | TotVectorTP==mappingTP(2)) & (TotNearestTarget==3 | TotNearestTarget==1),find(TimeVector==-200):find(TimeVector==50),2),1),'k:','LineWidth',2.5);
% xline(0,'k-.','Force onset','LineWidth',2.5);
xlabel('Time [ms]','FontSize',16); ylabel('Forward velocity [m/s]','FontSize',16); xlim([-200 50]);

subplot(1,2,2); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
plot(-200:50,mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(4) | TotVectorTP==mappingTP(2)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==50),2),1),'k','LineWidth',2.5);
plot(-200:50,mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(6) | TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(4) | TotVectorTP==mappingTP(2)) & (TotNearestTarget==3 | TotNearestTarget==1) & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-200):find(TimeVector==50),2),1),'k:','LineWidth',2.5);
% xline(0,'k-.','Force onset','LineWidth',2.5);
xlabel('Time [ms]','FontSize',16); ylabel('EMG activity [a.u]','FontSize',16); xlim([-200 50]);
%%
% figure for neuromatch 3 

close all; 
figure('units','normalized','outerposition',[0 0 1 1]); hold on;
subplot(2,2,1); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none'); set(gca,'XColor','none');
bar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1)],'ShowBaseLine','off','FaceColor',[1 0 0],'EdgeColor','none');  
xlim([0.3 2.7]); set(gca,'YColor','none');
er = errorbar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1)],[std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2)),std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2))]/sqrt(nSubjects),[std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2)),std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2))]/sqrt(nSubjects));
er.Color = [0 0 0]; er.LineWidth=3;

subplot(2,2,2); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');
bar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1)],'ShowBaseLine','off','FaceColor',[1 0 0],'EdgeColor','none');  
xlim([0.3 2.7]); set(gca,'YAxisLocation','right');ylabel('EMG activity [a.u.]');
er = errorbar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2),1)],[std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2)),std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2))]/sqrt(nSubjects),[std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2)),std(mean(TotEMG((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2))]/sqrt(nSubjects));
er.Color = [0 0 0]; er.LineWidth=3;

subplot(2,2,3); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');
bar([1,2],[mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1)],'ShowBaseLine','off','FaceColor',[1 0 0],'EdgeColor','none');  
xlim([0.3 2.7]); set(gca,'YColor','none');
er = errorbar([1 2],[mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1),mean(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2),1)],[std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2)),std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),1),2))]/sqrt(nSubjects),[std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2)),std(mean(TotEMG((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & (TotSubjects~=5 & TotSubjects~=10 & TotSubjects~=16 & TotSubjects~=18 & TotSubjects~=20),find(TimeVector==-150):find(TimeVector==0),2),2))]/sqrt(nSubjects));
er.Color = [0 0 0]; er.LineWidth=3;


%%



% The bottom part will contain information about the kinematics 
% position on the x-axis at the force onset
% forward speed at the force onset 

%TODO ajouter de la cinématique ci-dessous pour compléter la figure 
close all;
subplot(6,12,[52 53 54 64 65 66]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); title('un');
bar([1 2],[mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),1),mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==0),2),1)],'ShowBaseLine','off','EdgeColor','none','FaceColor',[1 0 0]); 
xlim([0.3 2.7]); xticks([1 2]); ylim([0 0.75]); xticklabels({'Center','Lateral'});
er = errorbar([1 2],[mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,TimeVector==0,2),1), mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,TimeVector==0,2),1)],[std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,TimeVector==0,2)), std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,TimeVector==0,2))]/sqrt(nSubjects),[std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,TimeVector==0,2)), std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,TimeVector==0,2))]/sqrt(nSubjects));
er.Color=[0 0 0]; er.LineWidth=3; ylabel({'Forward','speed [m/s]'});

subplot(6,12,[55 56 57 67 68 69]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18);set(gca,'Color','none'); set(gca,'YColor','none'); title('deux');
bar([1 2],[mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),1),mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),1)],'ShowBaseLine','off','EdgeColor','none','FaceColor',[1 0 0]); 
xlim([0.3 2.7]); xticks([1 2]); ylim([0 0.75]); xticklabels({'Center','Lateral'});
er = errorbar([1 2],[mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,TimeVector==0,2),1), mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,TimeVector==0,2),1)],[std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,TimeVector==0,2)), std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,TimeVector==0,2))]/sqrt(nSubjects),[std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,TimeVector==0,2)), std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,TimeVector==0,2))]/sqrt(nSubjects));
er.Color=[0 0 0]; er.LineWidth=3;
%%
% Let's do a better representation of the speed stuff
close all;
figure('Name','representation of the speed stuff','units','normalized'); hold on;
subplot(6,12,[52 53 54 64 65 66]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');
% for ii = 1 : nSubjects
%     plot([1,2], [mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2 & TotSubjects==ii,find(TimeVector==0),2),1),mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3 & TotSubjects==ii,find(TimeVector==0),2),1)],'.-','Color',[0.7 0.7 0.7],'LineWidth',2,'MarkerSize',25);
% end
plot([1,2], [mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),1),mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==0),2),1)],'k.-','LineWidth',3,'MarkerSize',35);
plot([1,1], [mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),1)+std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),0,1)/sqrt(nSubjects),mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),1)-std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),0,1)/sqrt(nSubjects)],'k','LineWidth',3);
plot([2,2], [mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==0),2),1)+std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==0),2),0,1)/sqrt(nSubjects),mean(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==0),2),1)-std(TotSpeed((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),2),0,1)/sqrt(nSubjects)],'k','LineWidth',3);
xlim([0.7 2.3]); ylim([0.35 0.45]);

subplot(6,12,[55 56 57 67 68 69]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YColor','none');
% for ii = 1 : nSubjects
%     plot([1,2], [mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2 & TotSubjects==ii,find(TimeVector==0),2),1),mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1 & TotSubjects==ii,find(TimeVector==0),2),1)],'.-','Color',[0.7 0.7 0.7],'LineWidth',2,'MarkerSize',25);
% end
plot([1,2], [mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),1),mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),1)],'k.-','LineWidth',3,'MarkerSize',35);
plot([1,1], [mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),1)+std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),0,1)/sqrt(nSubjects),mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),1)-std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),2),0,1)/sqrt(nSubjects)],'k','LineWidth',3);
plot([2,2], [mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),1)+std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),0,1)/sqrt(nSubjects),mean(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),1)-std(TotSpeed((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),2),0,1)/sqrt(nSubjects)],'k','LineWidth',3);
xlim([0.7 2.3]); ylim([0.35 0.45]);


%%


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

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(6,12,[55 56 57 67 68 69]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18);set(gca,'Color','none'); set(gca,'YColor','none');
bar([1,2],[mean(TotKine((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,find(TimeVector==0),1),1)-0.09,mean(TotKine((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,find(TimeVector==0),1),1)-0.09],'ShowBaseLine','off','EdgeColor','none','FaceColor',[1 0 0]);
xlim([0.3 2.7]); xticks([1 2]); ylim([-0.003 0.002]);
er = errorbar([1 2],[mean(TotKine((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,TimeVector==0,1),1), mean(TotKine((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,TimeVector==0,1),1)]-0.09,[std(TotKine((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,TimeVector==0,1)), std(TotKine((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,TimeVector==0,1))]/sqrt(nSubjects),[std(TotKine((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==2,TimeVector==0,1)), std(TotKine((TotVectorTP==mappingTP(8) | TotVectorTP==mappingTP(4)) & TotNearestTarget==3,TimeVector==0,1))]/sqrt(nSubjects));
er.Color=[0 0 0]; er.LineWidth=3;

subplot(6,12,[58 59 60 70 71 72]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18);set(gca,'Color','none'); set(gca,'YAxisLocation','right');
bar([1,2],[mean(TotKine((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,find(TimeVector==0),1),1)-0.09,mean(TotKine((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,find(TimeVector==0),1),1)-0.09],'ShowBaseLine','off','EdgeColor','none','FaceColor',[1 0 0]);
xlim([0.3 2.7]); xticks([1 2]); ylim([-0.003 0.002]);
er = errorbar([1 2],[mean(TotKine((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,TimeVector==0,1),1), mean(TotKine((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,TimeVector==0,1),1)]-0.09,[std(TotKine((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,TimeVector==0,1)), std(TotKine((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,TimeVector==0,1))]/sqrt(nSubjects),[std(TotKine((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==2,TimeVector==0,1)), std(TotKine((TotVectorTP==mappingTP(10) | TotVectorTP==mappingTP(6)) & TotNearestTarget==1,TimeVector==0,1))]/sqrt(nSubjects));
er.Color=[0 0 0]; er.LineWidth=3;
%%
% We'll run the stats here for the kinematics, it will be interesting as
% well

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





%% This section aims at giving an insight into the tradeoff between vigor and flexibility
close all;
% Computation of the speed at force onset 
SpeedZero = sqrt(TotSpeed(:,TimeVector==0,1).^2 + TotSpeed(:,TimeVector==0,2).^2);
SpeedZero = (SpeedZero - mean(SpeedZero)) / std(SpeedZero);

BinnedEMGValues = (BinnedEMGValues - mean(reshape(BinnedEMGValues,1,numel(BinnedEMGValues)),'omitnan'))/std(reshape(BinnedEMGValues,1,numel(BinnedEMGValues)),'omitnan');

figure('Name','Tradeoff vigor - flex','units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
for ii = 1 :nSubjects
    plot(mean(SpeedZero(TotSubjects==ii & TotNearestTarget==2)),mean(reshape(BinnedEMGValues([6,10],1,2,1,ii),1,numel(BinnedEMGValues([6,10],1,2,1,ii))),'omitnan'),'r.','MarkerSize',20);
    plot(mean(SpeedZero(TotSubjects==ii & (TotNearestTarget==1 | TotNearestTarget==3))),mean(reshape(BinnedEMGValues([6,10],1,[1,3],1,ii),1,numel(BinnedEMGValues([6,10],1,[1,3],1,ii))),'omitnan'),'b.','MarkerSize',20);
end

subplot(2,2,2); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
for ii = 1 :nSubjects
    plot(mean(SpeedZero(TotSubjects==ii & TotNearestTarget==2)),mean(reshape(BinnedEMGValues([6,10],1,2,2,ii),1,numel(BinnedEMGValues([6,10],1,2,2,ii))),'omitnan'),'r.','MarkerSize',20);
    plot(mean(SpeedZero(TotSubjects==ii & (TotNearestTarget==1 | TotNearestTarget==3))),mean(reshape(BinnedEMGValues([6,10],1,[1,3],2,ii),1,numel(BinnedEMGValues([6,10],1,[1,3],2,ii))),'omitnan'),'b.','MarkerSize',20);
end

subplot(2,2,3); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
for ii = 1 :nSubjects
    plot(mean(SpeedZero(TotSubjects==ii & TotNearestTarget==2)),mean(reshape(BinnedEMGValues([4,8],1,2,1,ii),1,numel(BinnedEMGValues([4,8],1,2,1,ii))),'omitnan'),'r.','MarkerSize',20);
    plot(mean(SpeedZero(TotSubjects==ii & (TotNearestTarget==1 | TotNearestTarget==3))),mean(reshape(BinnedEMGValues([4,8],1,[1,3],1,ii),1,numel(BinnedEMGValues([4,8],1,[1,3],1,ii))),'omitnan'),'b.','MarkerSize',20);
end

subplot(2,2,4); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none');
for ii = 1 :nSubjects
    plot(mean(SpeedZero(TotSubjects==ii & TotNearestTarget==2)),mean(reshape(BinnedEMGValues([4,8],1,2,2,ii),1,numel(BinnedEMGValues([4,8],1,2,2,ii))),'omitnan'),'r.','MarkerSize',20);
    plot(mean(SpeedZero(TotSubjects==ii & (TotNearestTarget==1 | TotNearestTarget==3))),mean(reshape(BinnedEMGValues([4,8],1,[1,3],2,ii),1,numel(BinnedEMGValues([4,8],1,[1,3],2,ii))),'omitnan'),'b.','MarkerSize',20);
end


%% Investigation of the EMG activity at movement onset  
close all;
figure('Name','hello');
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

figure('Name','Modification post round 1','units','normalized','outerposition',[0 0 1 1]); hold on;
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