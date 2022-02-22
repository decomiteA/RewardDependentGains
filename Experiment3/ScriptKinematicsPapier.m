% The aim of this script is to perform kinematics analysis for the reward
% experiment 3

clc; close all;
%% 0. Loading the data

% replace the following line with the path & name of the .mat file 
DataManip3 = load('DSExpe3.mat');
DataManip3 = DataManip3.ds3;
nSubjects = length(fieldnames(DataManip3));
Subjects = cell(nSubjects,1);
TotSubjects = zeros(nSubjects*480,1);
sizeTimeMatrix = zeros(nSubjects,1);
for ii = 1 : nSubjects
    Subjects{ii} = strcat('S',num2str(ii));
    TotSubjects((ii-1)*480+1:ii*480) = ii*ones(480,1);
    sizeTimeMatrix(ii) = size(DataManip3.(Subjects{ii}).matrix_timing,2);
end
extremeMin = DataManip3.(Subjects{1}).time_matrix(1);
extremeMax = DataManip3.(Subjects{1}).time_matrix(end);
sizeTimeMatrixTot = max(sizeTimeMatrix);
for ii = 2 : nSubjects
    extremeMin = max(extremeMin, DataManip3.(Subjects{ii}).time_matrix(1));
    if extremeMin==0
        disp(ii)
    end
    extremeMax = min(extremeMax, DataManip3.(Subjects{ii}).time_matrix(end));
end
TimeVector = extremeMin : extremeMax;
TotVectorTP = zeros(nSubjects*480,1);

TotReached = zeros(nSubjects*480,1);
TotSuccess = zeros(nSubjects*480,1);
TotKine = zeros(480*nSubjects,length(TimeVector),2);
TotTimeMat = zeros(480*nSubjects,sizeTimeMatrixTot);
idx_reached = zeros(480*nSubjects,1);
idx_success = zeros(480*nSubjects,1);
TotEMG = zeros(480*nSubjects,length(TimeVector),2);
for ii = 1 : nSubjects
    [TotReached((ii-1)*480+1:ii*480), idx_reached((ii-1)*480+1:ii*480)] = BooleanReachedM2Bis(DataManip3.(Subjects{ii}));
    TotVectorTP((ii-1)*480+1:ii*480) = DataManip3.(Subjects{ii}).vector_TP;
    [TotSuccess((ii-1)*480+1:ii*480),idx_success((ii-1)*480+1:ii*480)] = BooleanSuccessM2Bis(DataManip3.(Subjects{ii}));
    TotKine((ii-1)*480+1:ii*480,:,:) = DataManip3.(Subjects{ii}).matrix_kinematics(:,find(DataManip3.(Subjects{ii}).time_matrix==TimeVector(1)):find(DataManip3.(Subjects{ii}).time_matrix==TimeVector(end)),:);
    TotTimeMat((ii-1)*480+1:ii*480,1:sizeTimeMatrix(ii)) = DataManip3.(Subjects{ii}).matrix_timing;
    TotEMG((ii-1)*480+1:ii*480,:,:) = DataManip3.(Subjects{ii}).matrix_EMG(:,find(DataManip3.(Subjects{ii}).time_matrix==TimeVector(1)):find(DataManip3.(Subjects{ii}).time_matrix==TimeVector(end)),:);
end

MeanCoco = zeros(size(TotEMG,1),2);
TotSubRandom = zeros(size(TotEMG,1),1);
for ii = 1 : size(MeanCoco,1)
    MeanCoco(ii,:) = mean(TotEMG(ii,find(TimeVector==-500):find(TimeVector==-300),:),2);
    if (TotSubjects(ii)~=3 && TotSubjects(ii)~=6 && TotSubjects(ii)~=8 && TotSubjects(ii)~=11)
        TotSubRandom(ii) = TotSubjects(ii);
    end
end
mappingTP = [1 2 4 5 7 8 10 11 13 14];

%% Kinematics analysis
mappingTP = [1 2 4 5 7 8 10 11 13 14];
TotKineSpeed = zeros(size(TotKine,1),size(TotKine,2),2);
TotKineSpeed(:,3:end-2,:) = (-TotKine(:,5:end,:)+8*TotKine(:,4:end-1,:)-8*TotKine(:,2:end-3,:)+TotKine(:,1:end-4,:))/(12*0.001);

MaxSpeedX = PeakVelocity(TotKineSpeed(:,find(TimeVector==-100):find(TimeVector==200),1),-100:200);
MaxSpeedY = PeakVelocity(TotKineSpeed(:,find(TimeVector==-100):find(TimeVector==200),2),-100:200);

% Determination of the nearest target & initial reaching angle

TotNearTarget = zeros(nSubjects*480,1);
pos_targets = [0.03 0.28; 0.09 0.28; 0.15 0.28];
idxend = zeros(nSubjects*480,1);
for ii = 1 : nSubjects
    vectrick = DataManip3.(Subjects{ii}).time_matrix;
    for jj = 1 : 480
        if (TotReached((ii-1)*480+jj)==1)
            idxend((ii-1)*480+jj) = find(DataManip3.(Subjects{ii}).time_matrix==(min(abs(floor(DataManip3.(Subjects{ii}).matrix_timing(jj,idx_reached((ii-1)*480+jj))*1000)-floor(DataManip3.(Subjects{ii}).matrix_timing(jj,idx_reached((ii-1)*480+jj)-3)*1000)),floor(vectrick(end)))));
            TotNearTarget((ii-1)*480+jj) = NearestTarget([DataManip3.(Subjects{ii}).matrix_kinematics(jj,idxend((ii-1)*480+jj),1), DataManip3.(Subjects{ii}).matrix_kinematics(jj,idxend((ii-1)*480+jj),2)],pos_targets);
        end
    end
end
TotAngle = zeros(nSubjects*480,1);
for ii = 1 : nSubjects
    TotAngle((ii-1)*480+1:ii*480) = initialAngle(TotKine((ii-1)*480+1:ii*480,:,:),TimeVector);
end

TotSwitch = ~(TotNearTarget == 2 | TotNearTarget == 0);

ProportionTarget = zeros(10,3,nSubjects);
for ii = 1 : nSubjects
    for jj = 1 :10
        ProportionTarget(jj,:,ii) = [length(find(TotReached(TotSubjects==ii & TotVectorTP==mappingTP(jj) & TotNearTarget==1))) length(find(TotReached(TotSubjects==ii & TotVectorTP==mappingTP(jj) & TotNearTarget==2))) length(find(TotReached(TotSubjects==ii & TotVectorTP==mappingTP(jj) & TotNearTarget==3)))]/(length(find(TotReached(TotSubjects==ii & TotVectorTP==mappingTP(jj))==1)));
    end
end
%%
% Figure on switch proportion
id_random = zeros(nSubjects,1);
for ii = 1 : length(id_random)
    id_random(ii) = (ProportionTarget(2,2,ii)<0.25);
end
figure('Name','Investigation of the switch proportions','units','normalized','outerposition',[0 0 1 1]); hold on;

subplot(2,2,1); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none'); title('No background force - Left target');
for ii = 1 : nSubjects
    if id_random(ii)==0
        plot([1 2 3 4 5],[ProportionTarget(5,1,ii) ProportionTarget(9,1,ii) ProportionTarget(1,1,ii) ProportionTarget(7,1,ii) ProportionTarget(3,1,ii)],'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
    end
end
plot([1 2 3 4 5],[mean(ProportionTarget(5,1,id_random==0),3) mean(ProportionTarget(9,1,id_random==0),3) mean(ProportionTarget(1,1,id_random==0),3) mean(ProportionTarget(7,1,id_random==0),3) mean(ProportionTarget(3,1,id_random==0),3)],'k.-','MarkerSize',40,'LineWidth',3);
xlabel('Force [N]'); xlim([0.7 5.3]); ylim([-0.1 1.1]); yticks([0 0.25 0.5 0.75 1]); yticklabels({'0','25 %','50 %','75 %','100 %'}); ylabel('Reach proportion');
xticks([1 2 3 4 5]); xticklabels({'-6N','-3N','0N','3N','6N'});

subplot(2,2,2); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none'); title('No background force - Right target');
for ii = 1 : nSubjects
    if id_random(ii)==0
        plot([1 2 3 4 5],[ProportionTarget(5,3,ii) ProportionTarget(9,3,ii) ProportionTarget(1,3,ii) ProportionTarget(7,3,ii) ProportionTarget(3,3,ii)],'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
    end
end
plot([1 2 3 4 5],[mean(ProportionTarget(5,3,id_random==0),3) mean(ProportionTarget(9,3,id_random==0),3) mean(ProportionTarget(1,3,id_random==0),3) mean(ProportionTarget(7,3,id_random==0),3) mean(ProportionTarget(3,3,id_random==0),3)],'k.-','MarkerSize',40,'LineWidth',3);
xlabel('Force [N]'); xlim([0.7 5.3]); ylim([-0.1 1.1]); yticks([0 0.25 0.5 0.75 1]); yticklabels({'0','25 %','50 %','75 %','100 %'}); ylabel('Reach proportion');
xticks([1 2 3 4 5]); xticklabels({'-6N','-3N','0N','3N','6N'});

subplot(2,2,3); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none'); title('Background force - Left target');
for ii = 1 : nSubjects
    if id_random(ii)==0
        plot([1 2 3 4 5],[ProportionTarget(6,1,ii) ProportionTarget(10,1,ii) ProportionTarget(2,1,ii) ProportionTarget(8,1,ii) ProportionTarget(4,1,ii)],'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
    end
end
plot([1 2 3 4 5],[mean(ProportionTarget(6,1,id_random==0),3) mean(ProportionTarget(10,1,id_random==0),3) mean(ProportionTarget(2,1,id_random==0),3) mean(ProportionTarget(8,1,id_random==0),3) mean(ProportionTarget(4,1,id_random==0),3)],'k.-','MarkerSize',40,'LineWidth',3);
xlabel('Force [N]'); xlim([0.7 5.3]); ylim([-0.1 1.1]); yticks([0 0.25 0.5 0.75 1]); yticklabels({'0','25 %','50 %','75 %','100 %'}); ylabel('Reach proportion');
xticks([1 2 3 4 5]); xticklabels({'-6N','-3N','0N','3N','6N'});

subplot(2,2,4); hold on; set(gca,'LineWidth',2); set(gca,'FontSize',14); set(gca,'Color','none'); title('Background force - Right target');
for ii = 1 : nSubjects
    if id_random(ii)==0
        plot([1 2 3 4 5],[ProportionTarget(6,3,ii) ProportionTarget(10,3,ii) ProportionTarget(2,3,ii) ProportionTarget(8,3,ii) ProportionTarget(4,3,ii)],'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
    end
end
plot([1 2 3 4 5],[mean(ProportionTarget(6,3,id_random==0),3) mean(ProportionTarget(10,3,id_random==0),3) mean(ProportionTarget(2,3,id_random==0),3) mean(ProportionTarget(8,3,id_random==0),3) mean(ProportionTarget(4,3,id_random==0),3)],'k.-','MarkerSize',40,'LineWidth',3);
xlabel('Force [N]'); xlim([0.7 5.3]); ylim([-0.1 1.1]); yticks([0 0.25 0.5 0.75 1]); yticklabels({'0','25 %','50 %','75 %','100 %'}); ylabel('Reach proportion');
xticks([1 2 3 4 5]); xticklabels({'-6N','-3N','0N','3N','6N'});

[p56,tbl56,stats56] = anova1([reshape(ProportionTarget(5,1,id_random==0),nSubjects-2,1),reshape(ProportionTarget(6,1,id_random==0),nSubjects-2,1)],[],'off');
[p78,tbl78,stats78] = anova1([reshape(ProportionTarget(7,3,id_random==0),nSubjects-2,1),reshape(ProportionTarget(8,3,id_random==0),nSubjects-2,1)],[],'off');
[p34,tbl34,stats34] = anova1([reshape(ProportionTarget(3,3,id_random==0),nSubjects-2,1),reshape(ProportionTarget(4,3,id_random==0),nSubjects-2,1)],[],'off');
[p91,tbl91,stats91] = anova1([reshape(ProportionTarget(9,1,id_random==0),nSubjects-2,1),reshape(ProportionTarget(10,1,id_random==0),nSubjects-2,1)],[],'off');

p56_rsl = ranksum(reshape(ProportionTarget(5,1,id_random==0),nSubjects-2,1),reshape(ProportionTarget(6,1,id_random==0),nSubjects-2,1));
p78_rsl = ranksum(reshape(ProportionTarget(7,1,id_random==0),nSubjects-2,1),reshape(ProportionTarget(8,1,id_random==0),nSubjects-2,1));
p12_rsl = ranksum(reshape(ProportionTarget(1,1,id_random==0),nSubjects-2,1),reshape(ProportionTarget(2,1,id_random==0),nSubjects-2,1));
p34_rsl = ranksum(reshape(ProportionTarget(3,1,id_random==0),nSubjects-2,1),reshape(ProportionTarget(4,1,id_random==0),nSubjects-2,1));
p91_rsl = ranksum(reshape(ProportionTarget(9,1,id_random==0),nSubjects-2,1),reshape(ProportionTarget(10,1,id_random==0),nSubjects-2,1));
p56_rsr = ranksum(reshape(ProportionTarget(5,3,id_random==0),nSubjects-2,1),reshape(ProportionTarget(6,3,id_random==0),nSubjects-2,1));
p78_rsr = ranksum(reshape(ProportionTarget(7,3,id_random==0),nSubjects-2,1),reshape(ProportionTarget(8,3,id_random==0),nSubjects-2,1));
p12_rsr = ranksum(reshape(ProportionTarget(1,3,id_random==0),nSubjects-2,1),reshape(ProportionTarget(2,3,id_random==0),nSubjects-2,1));
p34_rsr = ranksum(reshape(ProportionTarget(3,3,id_random==0),nSubjects-2,1),reshape(ProportionTarget(4,3,id_random==0),nSubjects-2,1));
p91_rsr = ranksum(reshape(ProportionTarget(9,3,id_random==0),nSubjects-2,1),reshape(ProportionTarget(10,3,id_random==0),nSubjects-2,1));

%% Multilinear logistic regression

TotForce = zeros(nSubjects*480,1); TotCdt = mod(TotVectorTP,3);
TotSubRandom = zeros(nSubjects*480,1);
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
    if (TotSubjects(ii)~=3 & TotSubjects(ii)~=8)
        TotSubRandom(ii) = TotSubjects(ii);
    end
end
for ii = 1 : length(TotNearTarget)
    if TotNearTarget(ii)==2
        TotNearTarget(ii) = 4;
    end
end

%% generalized linear model with mixed effects 


[B,dev,stats] = mnrfit([TotForce(TotSubRandom~=0) TotCdt(TotSubRandom~=0),TotSubRandom(TotSubRandom~=0)],categorical(TotNearTarget(TotSubRandom~=0)),'Interactions','on');


tableglme = table((TotNearTarget(TotSubRandom~=0)),TotForce(TotSubRandom~=0),TotCdt(TotSubRandom~=0),TotSubRandom(TotSubRandom~=0),'VariableNames',{'target','force','background','subject'});
mrnglme = fitglme(tableglme,'target~force+background+(1|subject)','Distribution','Binomial','Link','logit','BinomialSize',length(TotNearTarget(TotSubRandom~=0)));


%%
TotNearTarget(TotNearTarget==2)=4;


nbtst_vec = 100:100:100; pvalbtst = zeros('like',nbtst_vec);
for jj = 1 : length(nbtst_vec)
nbtst = nbtst_vec(jj);
Bleft = zeros(nbtst,1);
Bright = zeros(nbtst,1);
idxtot = 1:length(TotForce);
Blefttot = zeros(nbtst,1);
Brighttot = zeros(nbtst,1);
for ii = 1 : nbtst
   idxleft = datasample(idxtot(TotSubRandom~=0 & TotNearTarget~=3 & TotNearTarget~=0),length(idxtot(TotSubRandom~=0 & TotNearTarget~=3 & TotNearTarget~=0)),'Replace',true);
   [Bleft_int,~,~] = mnrfit([TotForce(idxleft) TotCdt(idxleft) TotSubRandom(idxleft)],categorical(TotNearTarget(idxleft)));
   idxright = datasample(idxtot(TotSubRandom~=0 & TotNearTarget~=1 & TotNearTarget~=0),length(idxtot(TotSubRandom~=0 & TotNearTarget~=1 & TotNearTarget~=0)),'Replace',true);
   [Bright_int,~,~] = mnrfit([TotForce(idxright) TotCdt(idxright) TotSubRandom(idxright)],categorical(TotNearTarget(idxright)));
    Bleft(ii) = Bleft_int(3);
    Bright(ii) = Bright_int(3);
    idxtotm = datasample(idxtot(TotSubRandom~=0),length(idxtot(TotSubRandom~=0)),'Replace',true);
    [Btot_int,~,~] = mnrfit([TotForce(idxtotm) TotCdt(idxtotm) TotSubRandom(idxtotm)], categorical(TotNearTarget(idxtotm)));
    Blefttot(ii) = Btot_int(3,2);
    Brighttot(ii) = Btot_int(3,3);
end
[pvalbtst(jj),~] = signrank(Bright,Bleft);
end
for ii = 1 : length(TotNearTarget)
    if TotNearTarget(ii)==4
        TotNearTarget(ii) = 2;
    end
end

%%Is there an effect of the background load on the forward velocity
table_back_speed = table(TotKineSpeed(:,find(TimeVector==0),2),TotCdt,TotSubjects,'VariableNames',{'FwdSpeed','Cdt','Subjects'});
lme_back_speed = fitlme(table_back_speed,'FwdSpeed~Cdt+(1|Subjects)');

%%



% Investigation of the forward peak velocity

TotPeakVelocity = zeros(length(TotVectorTP),1);
for ii = 1 : length(TotPeakVelocity)
    TotPeakVelocity(ii) = max(TotKineSpeed(ii,find(TimeVector==-200):find(TimeVector==200),2));
end


tablelme1 = table(TotPeakVelocity((TotVectorTP==7 | TotVectorTP==8) & TotNearTarget==2), TotVectorTP((TotVectorTP==7 | TotVectorTP==8) & TotNearTarget==2), TotSubjects((TotVectorTP==7 | TotVectorTP==8) & TotNearTarget==2),'VariableNames',{'PeakSpeed','TP','Subjects'});
tablelme2 = table(TotPeakVelocity((TotVectorTP==13 | TotVectorTP==14) & TotNearTarget==2), TotVectorTP((TotVectorTP==13 | TotVectorTP==14) & TotNearTarget==2), TotSubjects((TotVectorTP==13 | TotVectorTP==14) & TotNearTarget==2),'VariableNames',{'PeakSpeed','TP','Subjects'});
tablelme3 = table(TotPeakVelocity((TotVectorTP==1 | TotVectorTP==2) & TotNearTarget==2), TotVectorTP((TotVectorTP==1 | TotVectorTP==2) & TotNearTarget==2), TotSubjects((TotVectorTP==1 | TotVectorTP==2) & TotNearTarget==2),'VariableNames',{'PeakSpeed','TP','Subjects'});
tablelme4 = table(TotPeakVelocity((TotVectorTP==10 | TotVectorTP==11) & TotNearTarget==2), TotVectorTP((TotVectorTP==10 | TotVectorTP==11) & TotNearTarget==2), TotSubjects((TotVectorTP==10 | TotVectorTP==11) & TotNearTarget==2),'VariableNames',{'PeakSpeed','TP','Subjects'});
tablelme5 = table(TotPeakVelocity((TotVectorTP==4 | TotVectorTP==5) & TotNearTarget==2), TotVectorTP((TotVectorTP==4 | TotVectorTP==5) & TotNearTarget==2), TotSubjects((TotVectorTP==4 | TotVectorTP==5) & TotNearTarget==2),'VariableNames',{'PeakSpeed','TP','Subjects'});

lme1 = fitlme(tablelme1,'PeakSpeed~TP+(1|Subjects)');
lme2 = fitlme(tablelme2,'PeakSpeed~TP+(1|Subjects)');
lme3 = fitlme(tablelme3,'PeakSpeed~TP+(1|Subjects)');
lme4 = fitlme(tablelme4,'PeakSpeed~TP+(1|Subjects)');
lme5 = fitlme(tablelme5,'PeakSpeed~TP+(1|Subjects)');

MatrixNoBck = zeros(nSubjects,1); matind = zeros(nSubjects,1);
MatrixBck = zeros(nSubjects,1);
for ii = 1 : nSubjects
   MatrixNoBck(ii) = mean(TotPeakVelocity(TotCdt==1 & TotNearTarget==2 & TotSubjects==ii));
   MatrixBck(ii) = mean(TotPeakVelocity(TotCdt==2 & TotNearTarget==2 & TotSubjects==ii));
   matind(ii) = mean(TotPeakVelocity(TotNearTarget==2 & TotSubjects==ii));
end
MatrixNoBck = MatrixNoBck-matind;
MatrixBck = MatrixBck-matind;


%% Paper Figure 
close all; 
figure('Name','Figure experiment 3','units','normalized','outerposition',[0 0 1 1]); hold on;

subplot(6,12,[1 2 3 4 13 14 15 16]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');
for ii = 1 : nSubjects
    if id_random(ii)==0
        plot([1 2 3 4 5],[ProportionTarget(5,1,ii) ProportionTarget(9,1,ii) ProportionTarget(1,1,ii) ProportionTarget(7,1,ii) ProportionTarget(3,1,ii)],'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
    end
end
plot([1 2 3 4 5],[mean(ProportionTarget(5,1,id_random==0),3) mean(ProportionTarget(9,1,id_random==0),3) mean(ProportionTarget(1,1,id_random==0),3) mean(ProportionTarget(7,1,id_random==0),3) mean(ProportionTarget(3,1,id_random==0),3)],'k.-','MarkerSize',40,'LineWidth',3);
xlabel('Force [N]'); xlim([0.7 5.3]); ylim([-0.1 1.1]); yticks([0 0.25 0.5 0.75 1]); yticklabels({'0','25 %','50 %','75 %','100 %'}); ylabel('Reach proportion');
xticks([1 2 3 4 5]); xticklabels({'-6N','-3N','0N','3N','6N'});

subplot(6,12,[5 6 7 8 17 18 19 20]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YColor','none');set(gca,'XColor','none');
for ii = 1 : nSubjects
    if id_random(ii)==0
        plot([1 2 3 4 5],[ProportionTarget(5,3,ii) ProportionTarget(9,3,ii) ProportionTarget(1,3,ii) ProportionTarget(7,3,ii) ProportionTarget(3,3,ii)],'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
    end
end
plot([1 2 3 4 5],[mean(ProportionTarget(5,3,id_random==0),3) mean(ProportionTarget(9,3,id_random==0),3) mean(ProportionTarget(1,3,id_random==0),3) mean(ProportionTarget(7,3,id_random==0),3) mean(ProportionTarget(3,3,id_random==0),3)],'k.-','MarkerSize',40,'LineWidth',3);
xlabel('Force [N]'); xlim([0.7 5.3]); ylim([-0.1 1.1]); yticks([0 0.25 0.5 0.75 1]); yticklabels({'0','25 %','50 %','75 %','100 %'}); ylabel('Reach proportion');
xticks([1 2 3 4 5]); xticklabels({'-6N','-3N','0N','3N','6N'});

subplot(6,12,[25 26 27 28 37 38 39 40]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'Xcolor','none');
for ii = 1 : nSubjects
    if id_random(ii)==0
        plot([1 2 3 4 5],[ProportionTarget(6,1,ii) ProportionTarget(10,1,ii) ProportionTarget(2,1,ii) ProportionTarget(8,1,ii) ProportionTarget(4,1,ii)],'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
    end
end
plot([1 2 3 4 5],[mean(ProportionTarget(6,1,id_random==0),3) mean(ProportionTarget(10,1,id_random==0),3) mean(ProportionTarget(2,1,id_random==0),3) mean(ProportionTarget(8,1,id_random==0),3) mean(ProportionTarget(4,1,id_random==0),3)],'k.-','MarkerSize',40,'LineWidth',3);
xlabel('Force [N]'); xlim([0.7 5.3]); ylim([-0.1 1.1]); yticks([0 0.25 0.5 0.75 1]); yticklabels({'0','25 %','50 %','75 %','100 %'}); ylabel('Reach proportion');
xticks([1 2 3 4 5]); xticklabels({'-6N','-3N','0N','3N','6N'});


subplot(6,12,[29 30 31 32 41 42 43 44]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YColor','none'); set(gca,'XColor','none');
for ii = 1 : nSubjects
    if id_random(ii)==0
        plot([1 2 3 4 5],[ProportionTarget(6,3,ii) ProportionTarget(10,3,ii) ProportionTarget(2,3,ii) ProportionTarget(8,3,ii) ProportionTarget(4,3,ii)],'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
    end
end
plot([1 2 3 4 5],[mean(ProportionTarget(6,3,id_random==0),3) mean(ProportionTarget(10,3,id_random==0),3) mean(ProportionTarget(2,3,id_random==0),3) mean(ProportionTarget(8,3,id_random==0),3) mean(ProportionTarget(4,3,id_random==0),3)],'k.-','MarkerSize',40,'LineWidth',3);
xlabel('Force [N]'); xlim([0.7 5.3]); ylim([-0.1 1.1]); yticks([0 0.25 0.5 0.75 1]); yticklabels({'0','25 %','50 %','75 %','100 %'}); ylabel('Reach proportion');
xticks([1 2 3 4 5]); xticklabels({'-6N','-3N','0N','3N','6N'});

subplot(6,12,[49 50 51 52 61 62 63 64]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); 
b1 = bar([1 2 4 5 7 8 10 11 13 14],[mean(ProportionTarget(5,1,id_random==0)) mean(ProportionTarget(6,1,id_random==0)) mean(ProportionTarget(9,1,id_random==0)) mean(ProportionTarget(10,1,id_random==0)) mean(ProportionTarget(1,1,id_random==0)) mean(ProportionTarget(2,1,id_random==0)) mean(ProportionTarget(7,1,id_random==0)) mean(ProportionTarget(8,1,id_random==0)) mean(ProportionTarget(3,1,id_random==0)) mean(ProportionTarget(4,1,id_random==0))],'FaceColor',[1 0 0],'EdgeColor','none','ShowBaseLine','off'); 
b1.FaceColor='flat';
b1.CData(2,:) = [0 0 1]; b1.CData(4,:) = [0 0 1]; b1.CData(6,:) = [0 0 1]; b1.CData(8,:) = [0 0 1]; b1.CData(10,:) = [0 0 1];
er1 = errorbar([1 2],[mean(ProportionTarget(5,1,id_random==0)), mean(ProportionTarget(6,1,id_random==0))],[std(ProportionTarget(5,1,id_random==0)), std(ProportionTarget(6,1,id_random==0))]/sqrt(length(find(id_random==0))),[std(ProportionTarget(5,1,id_random==0)), std(ProportionTarget(6,1,id_random==0))]/sqrt(length(find(id_random==0))));
er1.Color = [0 0 0]; er1.LineWidth=3;

er2 = errorbar([4 5],[mean(ProportionTarget(9,1,id_random==0)), mean(ProportionTarget(10,1,id_random==0))],[std(ProportionTarget(9,1,id_random==0)), std(ProportionTarget(10,1,id_random==0))]/sqrt(length(find(id_random==0))),[std(ProportionTarget(9,1,id_random==0)), std(ProportionTarget(10,1,id_random==0))]/sqrt(length(find(id_random==0))));
er2.Color = [0 0 0]; er2.LineWidth=3;

er3 = errorbar([7 8],[mean(ProportionTarget(1,1,id_random==0)), mean(ProportionTarget(2,1,id_random==0))],[std(ProportionTarget(1,1,id_random==0)), std(ProportionTarget(2,1,id_random==0))]/sqrt(length(find(id_random==0))),[std(ProportionTarget(1,1,id_random==0)), std(ProportionTarget(2,1,id_random==0))]/sqrt(length(find(id_random==0))));
er3.Color = [0 0 0]; er3.LineWidth=3;

er4 = errorbar([10 11],[mean(ProportionTarget(7,1,id_random==0)), mean(ProportionTarget(8,1,id_random==0))],[std(ProportionTarget(7,1,id_random==0)), std(ProportionTarget(8,1,id_random==0))]/sqrt(length(find(id_random==0))),[std(ProportionTarget(7,1,id_random==0)), std(ProportionTarget(8,1,id_random==0))]/sqrt(length(find(id_random==0))));
er4.Color = [0 0 0]; er4.LineWidth=3;

er5 = errorbar([13 14],[mean(ProportionTarget(3,1,id_random==0)), mean(ProportionTarget(4,1,id_random==0))],[std(ProportionTarget(3,1,id_random==0)), std(ProportionTarget(4,1,id_random==0))]/sqrt(length(find(id_random==0))),[std(ProportionTarget(3,1,id_random==0)), std(ProportionTarget(4,1,id_random==0))]/sqrt(length(find(id_random==0))));
er5.Color = [0 0 0]; er5.LineWidth=3;
xlabel('Force [N]'); xticklabels({'-6','-3','0','3','6'});xticks([1.5 4.5 7.5 10.5 13.5]);
yticks([0 0.3 0.6]); yticklabels({'0%','30%','60%'});

subplot(6,12,[53 54 55 56 65 66 67 68]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YColor','none');
b2 = bar([1 2 4 5 7 8 10 11 13 14],[mean(ProportionTarget(5,3,id_random==0)) mean(ProportionTarget(6,3,id_random==0)) mean(ProportionTarget(9,3,id_random==0)) mean(ProportionTarget(10,3,id_random==0)) mean(ProportionTarget(1,3,id_random==0)) mean(ProportionTarget(2,3,id_random==0)) mean(ProportionTarget(7,3,id_random==0)) mean(ProportionTarget(8,3,id_random==0)) mean(ProportionTarget(3,3,id_random==0)) mean(ProportionTarget(4,3,id_random==0))],'FaceColor',[1 0 0],'EdgeColor','none','ShowBaseLine','off'); 
b2.FaceColor='flat';
b2.CData(2,:) = [0 0 1]; b2.CData(4,:) = [0 0 1]; b2.CData(6,:) = [0 0 1]; b2.CData(8,:) = [0 0 1]; b2.CData(10,:) = [0 0 1];
er1 = errorbar([1 2],[mean(ProportionTarget(5,3,id_random==0)), mean(ProportionTarget(6,3,id_random==0))],[std(ProportionTarget(5,3,id_random==0)), std(ProportionTarget(6,3,id_random==0))]/sqrt(length(find(id_random==0))),[std(ProportionTarget(5,3,id_random==0)), std(ProportionTarget(6,3,id_random==0))]/sqrt(length(find(id_random==0))));
er1.Color = [0 0 0]; er1.LineWidth=3;
xticks([1.5 4.5 7.5 10.5 13.5]); xticklabels({'-6','-3','0','3','6'});
er2 = errorbar([4 5],[mean(ProportionTarget(9,3,id_random==0)), mean(ProportionTarget(10,3,id_random==0))],[std(ProportionTarget(9,3,id_random==0)), std(ProportionTarget(10,3,id_random==0))]/sqrt(length(find(id_random==0))),[std(ProportionTarget(9,3,id_random==0)), std(ProportionTarget(10,3,id_random==0))]/sqrt(length(find(id_random==0))));
er2.Color = [0 0 0]; er2.LineWidth=3;

er3 = errorbar([7 8],[mean(ProportionTarget(1,3,id_random==0)), mean(ProportionTarget(2,3,id_random==0))],[std(ProportionTarget(1,3,id_random==0)), std(ProportionTarget(2,3,id_random==0))]/sqrt(length(find(id_random==0))),[std(ProportionTarget(1,3,id_random==0)), std(ProportionTarget(2,3,id_random==0))]/sqrt(length(find(id_random==0))));
er3.Color = [0 0 0]; er3.LineWidth=3;

er4 = errorbar([10 11],[mean(ProportionTarget(7,3,id_random==0)), mean(ProportionTarget(8,3,id_random==0))],[std(ProportionTarget(7,3,id_random==0)), std(ProportionTarget(8,3,id_random==0))]/sqrt(length(find(id_random==0))),[std(ProportionTarget(7,3,id_random==0)), std(ProportionTarget(8,3,id_random==0))]/sqrt(length(find(id_random==0))));
er4.Color = [0 0 0]; er4.LineWidth=3;

er5 = errorbar([13 14],[mean(ProportionTarget(3,3,id_random==0)), mean(ProportionTarget(4,3,id_random==0))],[std(ProportionTarget(3,3,id_random==0)), std(ProportionTarget(4,3,id_random==0))]/sqrt(length(find(id_random==0))),[std(ProportionTarget(3,3,id_random==0)), std(ProportionTarget(4,3,id_random==0))]/sqrt(length(find(id_random==0))));
er5.Color = [0 0 0]; er5.LineWidth=3;



subplot(6,12,[9 10 11 12 21 22 23 24]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YAxisLocation','right');
xticks([1.5 4.5]); xticklabels({'PM','PD'}); ylabel('EMG activity [u.a.]');
b3 = bar([1 2 4 5],[mean(MeanCoco(TotCdt==1 & TotSubRandom~=0),1) mean(MeanCoco(TotCdt==2 & TotSubRandom~=0),1) mean(MeanCoco(TotCdt==1 & TotSubRandom~=0,2)) mean(MeanCoco(TotCdt==2 & TotSubRandom~=0,2))],'FaceColor',[1 0 0],'ShowBaseLine','off','EdgeColor','none');
b3.FaceColor = 'flat';
b3.CData(2,:)=[0 0 1]; b3.CData(4,:)=[0 0 1];
er1 = errorbar([1 2],[mean(MeanCoco(TotCdt==1 & TotSubRandom~=0,1)) mean(MeanCoco(TotCdt==2 & TotSubRandom~=0,1))],[mean(MeanCoco(TotCdt==1 & TotSubRandom~=0,1)) mean(MeanCoco(TotCdt==2 & TotSubRandom~=0,1))]/sqrt(15),[mean(MeanCoco(TotCdt==1 & TotSubRandom~=0,1)) mean(MeanCoco(TotCdt==2 & TotSubRandom~=0,1))]/sqrt(15));
er1.Color = [0 0 0]; er1.LineWidth=3;

er2 = errorbar([4 5],[mean(MeanCoco(TotCdt==1 & TotSubRandom~=0,2)) mean(MeanCoco(TotCdt==2 & TotSubRandom~=0,2))],[mean(MeanCoco(TotCdt==1 & TotSubRandom~=0,2)) mean(MeanCoco(TotCdt==2 & TotSubRandom~=0,2))]/sqrt(15),[mean(MeanCoco(TotCdt==1 & TotSubRandom~=0,2)) mean(MeanCoco(TotCdt==2 & TotSubRandom~=0,2))]/sqrt(15));
er2.Color = [0 0 0]; er2.LineWidth=3;

tableCoco_pec = table(MeanCoco((TotCdt==1 | TotCdt==2) & TotSubRandom~=0,1),TotCdt((TotCdt==1 | TotCdt==2) & TotSubRandom~=0), TotSubRandom((TotCdt==1 | TotCdt==2) & TotSubRandom~=0),'VariableNames',{'Coco','Condition','Subjects'});
tableCoco_del = table(MeanCoco((TotCdt==1 | TotCdt==2) & TotSubRandom~=0,2),TotCdt((TotCdt==1 | TotCdt==2) & TotSubRandom~=0), TotSubRandom((TotCdt==1 | TotCdt==2) & TotSubRandom~=0),'VariableNames',{'Coco','Condition','Subjects'});
lmecoco_pec = fitlme(tableCoco_pec,'Coco~Condition+(1|Subjects)');
lmecoco_del = fitlme(tableCoco_del,'Coco~Condition+(1|Subjects)');

subplot(6,12,[33 34 35 36 45 46 47 48]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YAxisLocation','right');
plot(-200:500, mean(TotKineSpeed(TotCdt==1 & TotNearTarget==2,find(TimeVector==-200):find(TimeVector==500),2),1),'r','LineWidth',3);
plot(-200:500, mean(TotKineSpeed(TotCdt==2 & TotNearTarget==2,find(TimeVector==-200):find(TimeVector==500),2),1),'b','LineWidth',3);
ylabel('Forward speed [m/s]'); ylim([0 1]); xlim([-220 520]); xticks([-200 0 500]); xticklabels({'-0.2','0','0.5'});

subplot(6,12,[57 58 59 60 69 70 71 72]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YAxisLocation','right');
b4 = bar([1 3],[mean(TotPeakVelocity(TotCdt==1 & TotNearTarget==2)) mean(TotPeakVelocity(TotCdt==2 & TotNearTarget==2))],'EdgeColor','none','FaceColor',[1 0 0],'ShowBaseLine','off','BarWidth',0.5); 
b4.FaceColor='flat'; b4.CData(2,:) = [0 0 1];
erfinal = errorbar([1 3],[mean(TotPeakVelocity(TotCdt==1 & TotNearTarget==2)) mean(TotPeakVelocity(TotCdt==2 & TotNearTarget==2))],[std(TotPeakVelocity(TotCdt==1 & TotNearTarget==2)) std(TotPeakVelocity(TotCdt==2 & TotNearTarget==2))]/sqrt(nSubjects),[std(TotPeakVelocity(TotCdt==1 & TotNearTarget==2)) std(TotPeakVelocity(TotCdt==2 & TotNearTarget==2))]/sqrt(nSubjects));
erfinal.LineWidth=3; erfinal.Color=[0 0 0]; xlim([0 4]); xticks([1 3]); xticklabels({'No background','Background'});
ylim([0 1]);

tablepeakvel = table(TotPeakVelocity((TotCdt==1 | TotCdt==2) & TotNearTarget==2),TotCdt((TotCdt==1 | TotCdt==2) & TotNearTarget==2),TotSubjects((TotCdt==1 | TotCdt==2) & TotNearTarget==2),'VariableNames',{'Speed','Condition','Subjects'});
lmepeakvel = fitlme(tablepeakvel,'Speed~Condition+(1|Subjects)');




%%
% Anova on muscle co-activation 

AnovaMatrix = zeros(7200,2);
TotSAnova = [1 2 4 5 7 9 10 12 13 14 15 16 17 18 19];
for ii=1:length(TotSAnova)
    AnovaMatrix((ii-1)*240+1:ii*240,:) = [MeanCoco(TotSubRandom==TotSAnova(ii) & mod(TotVectorTP,3)==1,1), MeanCoco(TotSubRandom==TotSAnova(ii) & mod(TotVectorTP,3)==1,2)];
    AnovaMatrix(3600+(ii-1)*240+1:3600+ii*240,:) = [MeanCoco(TotSubRandom==TotSAnova(ii) & mod(TotVectorTP,3)==2,1), MeanCoco(TotSubRandom==TotSAnova(ii) & mod(TotVectorTP,3)==2,2)];
end

[p_anova2,tbl_anova2,stats_anova2] = anova2(AnovaMatrix,3600);



% Comparison with a classical mixed model

table_lme_coco = table([MeanCoco(TotSubRandom~=0,1); MeanCoco(TotSubRandom~=0,2)],[TotSubRandom(TotSubRandom~=0);TotSubRandom(TotSubRandom~=0)],[ones(length(TotSubRandom(TotSubRandom~=0)),1);2*ones(length(TotSubRandom(TotSubRandom~=0)),1)],mod([TotVectorTP(TotSubRandom~=0);TotVectorTP(TotSubRandom~=0)],3),'VariableNames',{'EMG','Subjects','Muscle','Background'});
lme_coco = fitlme(table_lme_coco,'EMG~Muscle*Background+(1|Subjects)');


%%
% paired t-test individual tests
LeftNoB = zeros(5*length(id_random(id_random==0)),1);
LeftB = zeros(5*length(id_random(id_random==0)),1);
RightNoB = zeros(5*length(id_random(id_random==0)),1);
RightB = zeros(5*length(id_random(id_random==0)),1);
map1 = [5 9 1 7 3];
map2 = [6 10 2 8 4];
for ii = 1 : length(map1)
   LeftNoB((ii-1)*17+1:ii*17) = ProportionTarget(map1(ii),1,id_random==0);
   LeftB((ii-1)*17+1:ii*17) = ProportionTarget(map2(ii),1,id_random==0);
   RightNoB((ii-1)*17+1:ii*17) = ProportionTarget(map1(ii),3,id_random==0);
   RightB((ii-1)*17+1:ii*17) = ProportionTarget(map2(ii),3,id_random==0);
end

%% Investigation of the reaction time 


TotReactionTime = zeros(size(TotKine,1),1);
for ii = 1 : nSubjects
   TotReactionTime((ii-1)*480+1:480*ii) = DataManip3.(Subjects{ii}).matrix_timing(:,5) - DataManip3.(Subjects{ii}).matrix_timing(:,4); 
end

tableRT_cdt = table(TotReactionTime,TotCdt,TotSubjects,'VariableNames',{'RT','Condition','Subject'});
lmeRT_cdt = fitlme(tableRT_cdt,'RT~Condition+(1|Subject)');


tableRT_12 = table(TotReactionTime(TotCdt==1 & (TotNearTarget==1 | TotNearTarget==2)), TotNearTarget(TotCdt==1 & (TotNearTarget==1 | TotNearTarget==2)), TotSubjects(TotCdt==1 & (TotNearTarget==1 | TotNearTarget==2)),'VariableNames',{'RT','Target','Subject'});
lmeRT_12 = fitlme(tableRT_12,'RT~Target+(1|Subject)');

tableRT_23 = table(TotReactionTime(TotCdt==2 & (TotNearTarget==3 | TotNearTarget==2)), TotNearTarget(TotCdt==2 & (TotNearTarget==3 | TotNearTarget==2)), TotSubjects(TotCdt==2 & (TotNearTarget==3 | TotNearTarget==2)),'VariableNames',{'RT','Target','Subject'});
lmeRT_23 = fitlme(tableRT_23,'RT~Target+(1|Subject)');

%% Controlling for an effect of fatigue or habituation on the velocity
clc;
TotBlocks = repmat([1*ones(80,1);2*ones(80,1);3*ones(80,1);4*ones(80,1);5*ones(80,1);6*ones(80,1)],nSubjects,1);
table_fatigue_block = table(TotPeakVelocity(1:size(TotPeakVelocity,1)/2,:,:),TotBlocks(1:size(TotPeakVelocity,1)/2,:,:),TotSubjects(1:size(TotPeakVelocity,1)/2,:,:),'VariableNames',{'pv','block','subjects'});
lme_fatigue_block = fitlme(table_fatigue_block,'pv~block+(1|subjects)');