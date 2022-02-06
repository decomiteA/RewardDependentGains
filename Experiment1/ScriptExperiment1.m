% This scripts plots figures 2 and 3 and computes all the associated
% statistics of the following paper : 
% Reward-dependent selection of feedback gains impacts rapid motor decision
% De Comite A., Crevecoeur F., Lefèvre P. 
%
% Author : Antoine DE COMITE
% antoinedecomite@gmail.com
%
%
close all; clc; clear all;

%% 0. Data loading and pre-processing

M1Data = load('M1Data.mat');
nSubjects = length(fieldnames(M1Data));
Subjects = cell(nSubjects,1);
for ii = 1 : nSubjects
    Subjects{ii} = strcat('S',num2str(ii));
end
TotSubjects = zeros(nSubjects*432,1);
TotVectorTP = zeros(nSubjects*432,1);
for ii = 1 : nSubjects
    TotVectorTP((ii-1)*432+1:ii*432) = M1Data.(Subjects{ii}).vector_TP;
    TotSubjects((ii-1)*432+1:ii*432) = ii*ones(432,1);
end

extremeMin = M1Data.(Subjects{ii}).time_matrix(1);
extremeMax = M1Data.(Subjects{ii}).time_matrix(end);
for ii = 2 : nSubjects
    extremeMin = max(extremeMin,M1Data.(Subjects{ii}).time_matrix(1));
    extremeMax = min(extremeMax,M1Data.(Subjects{ii}).time_matrix(end));
end
TimeVector = extremeMin:extremeMax;                   % Time vector (common)
TotKine = zeros(nSubjects*432, length(TimeVector),2); % Raw kinematics data
TotEMG = zeros(nSubjects*432, length(TimeVector),2);  % Raw EMG data
for ii = 1 : nSubjects
    TotKine((ii-1)*432+1:ii*432,:,:) = M1Data.(Subjects{ii}).matrix_kinematics(:,find(M1Data.(Subjects{ii}).time_matrix==extremeMin):find(M1Data.(Subjects{ii}).time_matrix==extremeMax),:);
    TotEMG((ii-1)*432+1:ii*432,:,:) = M1Data.(Subjects{ii}).matrix_EMG(:,find(M1Data.(Subjects{ii}).time_matrix==extremeMin):find(M1Data.(Subjects{ii}).time_matrix==extremeMax),:);
end

% Subtracting kinematics baseline 
TotKineBaseline = zeros(nSubjects,length(TimeVector),3);
TotVectorTPDeviation = zeros(2016,1); TotSubjectsDeviation = zeros(2016,1);
for ii = 1 : nSubjects
    for jj = 1 : 3
        TotKineBaseline(ii,:,jj) = mean(TotKine(TotVectorTP==jj & TotSubjects==ii,:,1),1);
    end
    TotVectorTPDeviation((ii-1)*144+1:ii*144)=TotVectorTP(TotSubjects==ii & TotVectorTP~=1 & TotVectorTP~=2 & TotVectorTP~=3);
    TotSubjectsDeviation((ii-1)*144+1:ii*144)=TotSubjects(TotSubjects==ii & TotVectorTP~=1 & TotVectorTP~=2 & TotVectorTP~=3);
end

TotKineDeviation = zeros(2016,length(TimeVector));
count = 1; 
for ii = 1 : size(TotKine,1)
    if (TotVectorTP(ii)~=1 && TotVectorTP(ii)~=2 && TotVectorTP(ii)~=3)
        TotKineDeviation(count,:) = TotKine(ii,:,1) - TotKineBaseline(TotSubjects(ii),:,mappingMod(TotVectorTP(ii)));
        count = count + 1;
    end
end
TotKine(TotKine==0) = NaN;

% catching empty cells

for ii = 1 : nSubjects
    for jj = 1 :432
        if isempty(M1Data.(Subjects{ii}).array_timing{jj,8})
            M1Data.(Subjects{ii}).array_timing{jj,8} = 'buffer';
        end
    end
end


%% 1. Computation of the different metrics 

% a. Peak velocity

TotSpeed = zeros(size(TotKine,1),size(TotKine,2),2);
TotSpeed(:,3:end-2,:) = (-TotKine(:,5:end,:)+8*TotKine(:,4:end-1,:)-8*TotKine(:,2:end-3,:)+TotKine(:,1:end-4,:))/(12*0.001);
NormSpeed = sqrt(TotSpeed(:,:,1).^2 + TotSpeed(:,:,2).^2);

PeakVelocity = zeros(size(TotSpeed,1),1);
for ii = 1 : length(PeakVelocity)
    PeakVelocity(ii) = max(TotSpeed(ii,find(TimeVector==-200):find(TimeVector==200),2));
end

MatrixSpeedMeanTrial123 = zeros(nSubjects,3); MatrixSpeedMeanTrial456 = zeros(nSubjects,3); MatrixSpeedMeanTrial789 = zeros(nSubjects,3);
matpert1 = zeros(nSubjects,1); matpert4 = zeros(nSubjects,1); matpert7 = zeros(nSubjects,1);
for ii = 1 : size(MatrixSpeedMeanTrial123,1)
   MatrixSpeedMeanTrial123(ii,:) = [mean(PeakVelocity(TotSubjects==ii & TotVectorTP==1)), mean(PeakVelocity(TotSubjects==ii & TotVectorTP==2)), mean(PeakVelocity(TotSubjects==ii & TotVectorTP==3))];
   MatrixSpeedMeanTrial456(ii,:) = [mean(PeakVelocity(TotSubjects==ii & TotVectorTP==4)), mean(PeakVelocity(TotSubjects==ii & TotVectorTP==5)), mean(PeakVelocity(TotSubjects==ii & TotVectorTP==6))];
   MatrixSpeedMeanTrial789(ii,:) = [mean(PeakVelocity(TotSubjects==ii & TotVectorTP==7)), mean(PeakVelocity(TotSubjects==ii & TotVectorTP==8)), mean(PeakVelocity(TotSubjects==ii & TotVectorTP==9))];
   matpert1(ii) = mean(PeakVelocity(TotSubjects==ii & (TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3)));
   matpert4(ii) = mean(PeakVelocity(TotSubjects==ii & (TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6)));
   matpert7(ii) = mean(PeakVelocity(TotSubjects==ii & (TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9)));
end
std0 = std(reshape(MatrixSpeedMeanTrial123,1,numel(MatrixSpeedMeanTrial123)));
std1 = std(reshape(MatrixSpeedMeanTrial456,1,numel(MatrixSpeedMeanTrial456)));
std2 = std(reshape(MatrixSpeedMeanTrial456,1,numel(MatrixSpeedMeanTrial789)));
MatrixSpeedMeanTrial123 = MatrixSpeedMeanTrial123 - matpert1;
MatrixSpeedMeanTrial456 = MatrixSpeedMeanTrial456 - matpert4;
MatrixSpeedMeanTrial789 = MatrixSpeedMeanTrial789 - matpert7;
MatrixSpeedMeanTrial456z = MatrixSpeedMeanTrial456/std1;
MatrixSpeedMeanTrial789z = MatrixSpeedMeanTrial789/std2;
MatrixSpeedMeanTrial123z = MatrixSpeedMeanTrial123/std0;


% b. Hand deviation

MaxDeviation = max(abs(TotKineDeviation'));


MatrixHandTrial456 = zeros(nSubjects,3);MatrixHandTrial789 = zeros(nSubjects,3);
matpert4hand = zeros(nSubjects,1);matpert7hand = zeros(nSubjects,1);
for ii = 1 : size(MatrixHandTrial456,1)
    MatrixHandTrial456(ii,:) = [mean(MaxDeviation(TotVectorTPDeviation==4 & TotSubjectsDeviation==ii)), mean(MaxDeviation(TotVectorTPDeviation==5 & TotSubjectsDeviation==ii)), mean(MaxDeviation(TotVectorTPDeviation==6 & TotSubjectsDeviation==ii))];
    MatrixHandTrial789(ii,:) = [mean(MaxDeviation(TotVectorTPDeviation==7 & TotSubjectsDeviation==ii)), mean(MaxDeviation(TotVectorTPDeviation==8 & TotSubjectsDeviation==ii)), mean(MaxDeviation(TotVectorTPDeviation==9 & TotSubjectsDeviation==ii))];
    matpert4hand(ii) = mean(MaxDeviation(TotSubjectsDeviation==ii & (TotVectorTPDeviation==4 | TotVectorTPDeviation==5 | TotVectorTPDeviation==6)));
    matpert7hand(ii) = mean(MaxDeviation(TotSubjectsDeviation==ii & (TotVectorTPDeviation==7 | TotVectorTPDeviation==8 | TotVectorTPDeviation==9)));
end
std1 = std(reshape(MatrixHandTrial456,1,numel(MatrixHandTrial456)));
std2 = std(reshape(MatrixHandTrial789,1,numel(MatrixHandTrial789)));
MatrixHandTrial456 = MatrixHandTrial456 - matpert4hand;
MatrixHandTrial789 = MatrixHandTrial789 - matpert7hand;
MatrixHandTrial456z = MatrixHandTrial456/std1;
MatrixHandTrial789z = MatrixHandTrial789/std2;

% C. EMG activities - baseline
    
MatrixBinnedEMG = mean(TotEMG(:,find(TimeVector==0):find(TimeVector==200),:),2);
MatrixEMGTrial123Pec = zeros(nSubjects,3); MatrixEMGTrial123Del = zeros(nSubjects,3);
matmean1Pec = zeros(nSubjects,1); matmean1Del = zeros(nSubjects,1);
for ii = 1 : size(MatrixEMGTrial123Pec,1)
   MatrixEMGTrial123Pec(ii,:) = [mean(MatrixBinnedEMG(TotSubjects==ii & TotVectorTP==1,1)), mean(MatrixBinnedEMG(TotSubjects==ii & TotVectorTP==2,1)) mean(MatrixBinnedEMG(TotSubjects==ii & TotVectorTP==3,1))];  
   MatrixEMGTrial123Del(ii,:) = [mean(MatrixBinnedEMG(TotSubjects==ii & TotVectorTP==1,2)), mean(MatrixBinnedEMG(TotSubjects==ii & TotVectorTP==2,2)) mean(MatrixBinnedEMG(TotSubjects==ii & TotVectorTP==3,2))];  
   matmean1Pec(ii) = mean(MatrixBinnedEMG(TotSubjects==ii & (TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3),1));
   matmean1Del(ii) = mean(MatrixBinnedEMG(TotSubjects==ii & (TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3),2));
end
MatrixEMGTrial123Pecph = MatrixEMGTrial123Pec;
MatrixEMGTrial123Delph = MatrixEMGTrial123Del;
MatrixEMGTrial123Pec = MatrixEMGTrial123Pec - matmean1Pec;
MatrixEMGTrial123Del = MatrixEMGTrial123Del - matmean1Del;  
std1 = std(reshape(MatrixEMGTrial123Pecph,1,numel(MatrixEMGTrial123Pecph)));
std2 = std(reshape(MatrixEMGTrial123Delph,1,numel(MatrixEMGTrial123Delph)));
MatrixEMGTrial123Pecz = MatrixEMGTrial123Pec/std1;
MatrixEMGTrial123Delz = MatrixEMGTrial123Del/std2;

% d. EMG activities - feedback responses 

Bins = [-150   0;
         50  100;
         100 180];
    
MatrixBinnedREMG = zeros(size(TotEMG,1),3,2);
for ii = 1 : size(MatrixBinnedREMG,2)
   MatrixBinnedREMG(:,ii,:) = mean(TotEMG(:,find(TimeVector==Bins(ii,1)):find(TimeVector==Bins(ii,2)),:),2); 
end

% Here below - only the agonist responses 
MatrixEMGPec456 = zeros(nSubjects,2,3);MatrixEMGDel789 = zeros(nSubjects,2,3);
matmeanpec = zeros(nSubjects,2);matmeandel = zeros(nSubjects,2);
for ii = 1 : size(MatrixEMGPec456,1)
   MatrixEMGPec456(ii,1,:) = [mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==4,1,1)),mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==5,1,1)),mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==6,1,1))];
   MatrixEMGPec456(ii,2,:) = [mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==4,2,1)),mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==5,2,1)),mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==6,2,1))];
   MatrixEMGDel789(ii,1,:) = [mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==7,1,2)),mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==8,1,2)),mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==9,1,2))];
   MatrixEMGDel789(ii,2,:) = [mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==7,2,2)),mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==8,2,2)),mean(MatrixBinnedREMG(TotSubjects==ii & TotVectorTP==9,2,2))];
   matmeanpec(ii,1) = mean(MatrixBinnedREMG(TotSubjects==ii & (TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),1,1));
   matmeanpec(ii,2) = mean(MatrixBinnedREMG(TotSubjects==ii & (TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),2,1));
   matmeandel(ii,1) = mean(MatrixBinnedREMG(TotSubjects==ii & (TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9),1,2));
   matmeandel(ii,2) = mean(MatrixBinnedREMG(TotSubjects==ii & (TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9),2,2));
end
MatrixEMGPec456ph = MatrixEMGPec456;
MatrixEMGDel789ph = MatrixEMGDel789;
MatrixEMGPec456(:,1,:) = MatrixEMGPec456(:,1,:) - matmeanpec(:,1);
MatrixEMGPec456(:,2,:) = MatrixEMGPec456(:,2,:) - matmeanpec(:,2);
MatrixEMGDel789(:,1,:) = MatrixEMGDel789(:,1,:) - matmeandel(:,1);
MatrixEMGDel789(:,2,:) = MatrixEMGDel789(:,2,:) - matmeandel(:,2);
stda = std(reshape(MatrixEMGPec456ph(:,1,:),1,numel(MatrixEMGPec456ph(:,1,:))));
stdb = std(reshape(MatrixEMGPec456ph(:,2,:),1,numel(MatrixEMGPec456ph(:,2,:))));
stdc = std(reshape(MatrixEMGDel789ph(:,1,:),1,numel(MatrixEMGDel789ph(:,1,:))));
stdd = std(reshape(MatrixEMGDel789ph(:,2,:),1,numel(MatrixEMGDel789ph(:,2,:))));
MatrixEMGPec456LLz =  MatrixEMGPec456(:,1,:)/stda;
MatrixEMGPec456VOLz = MatrixEMGPec456(:,2,:)/stdb;
MatrixEMGDel789LLz =  MatrixEMGDel789(:,1,:)/stdc;
MatrixEMGDel789VOLz = MatrixEMGDel789(:,2,:)/stdd;



%% 2. Statistical tests for the different effects
    
% a. Peak velocity

tablepv123 = table(PeakVelocity(TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3),TotVectorTP(TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3),TotSubjects(TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3),'VariableNames',{'PeakVelocity','TP','Subjects'});
lmepv132 = fitlme(tablepv123,'PeakVelocity~TP+(1|Subjects)');
tablepv456 = table(PeakVelocity(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),TotVectorTP(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),TotSubjects(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),'VariableNames',{'PeakVelocity','TP','Subjects'});
lmepv456 = fitlme(tablepv456,'PeakVelocity~TP+(1|Subjects)');
tablepv789 = table(PeakVelocity(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9),TotVectorTP(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9),TotSubjects(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9),'VariableNames',{'PeakVelocity','TP','Subjects'});
lmepv789 = fitlme(tablepv789,'PeakVelocity~TP+(1|Subjects)');

% b. Hand deviation

tablehand_deviation456 = table(MaxDeviation(TotVectorTPDeviation==4 | TotVectorTPDeviation==5 | TotVectorTPDeviation==6)',TotVectorTPDeviation(TotVectorTPDeviation==4 | TotVectorTPDeviation==5 | TotVectorTPDeviation==6),TotSubjectsDeviation(TotVectorTPDeviation==4 | TotVectorTPDeviation==5 | TotVectorTPDeviation==6),'VariableNames',{'HandDeviation','Trial','Subjects'});
lmehand_deviation456 = fitlme(tablehand_deviation456,'HandDeviation~Trial+(1|Subjects)');
tablehand_deviation789 = table(MaxDeviation(TotVectorTPDeviation==7 | TotVectorTPDeviation==8 | TotVectorTPDeviation==9)',TotVectorTPDeviation(TotVectorTPDeviation==7 | TotVectorTPDeviation==8 | TotVectorTPDeviation==9),TotSubjectsDeviation(TotVectorTPDeviation==7 | TotVectorTPDeviation==8 | TotVectorTPDeviation==9),'VariableNames',{'HandDeviation','Trial','Subjects'});
lmehand_deviation789 = fitlme(tablehand_deviation789,'HandDeviation~Trial+(1|Subjects)');


% c. EMG activity - baseline

table123Pec = table(MatrixBinnedEMG(TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3,1), TotVectorTP(TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3),TotSubjects(TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3),'VariableNames',{'EMG','TP','Subjects'});
fitlme123Pec = fitlme(table123Pec,'EMG~TP+(1|Subjects)');
table123Del = table(MatrixBinnedEMG(TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3,2), TotVectorTP(TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3),TotSubjects(TotVectorTP==1 | TotVectorTP==2 | TotVectorTP==3),'VariableNames',{'EMG','TP','Subjects'});
fitlme123Del = fitlme(table123Del,'EMG~TP+(1|Subjects)');


% d. EMG activity - feedback responses

% Agonist responses 
table456PEC1 = table(MatrixBinnedREMG(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6,2,1), TotVectorTP(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),TotSubjects(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),'VariableNames',{'EMG','TP','Subjects'});
fitlme456PEC1 = fitlme(table456PEC1,'EMG~TP+(1|Subjects)');
table456PEC2 = table(MatrixBinnedREMG(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6,3,1), TotVectorTP(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),TotSubjects(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),'VariableNames',{'EMG','TP','Subjects'});
fitlme456PEC2 = fitlme(table456PEC2,'EMG~TP+(1|Subjects)');
table789Del1 = table(MatrixBinnedREMG(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9,2,2), TotVectorTP(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),TotSubjects(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),'VariableNames',{'EMG','TP','Subjects'});
fitlme789Del1 = fitlme(table789Del1,'EMG~TP+(1|Subjects)');
table789Del2 = table(MatrixBinnedREMG(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9,3,2), TotVectorTP(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),TotSubjects(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),'VariableNames',{'EMG','TP','Subjects'});
fitlme789Del2 = fitlme(table789Del2,'EMG~TP+(1|Subjects)');

% Antagonist responses
table456Del1 = table(MatrixBinnedREMG(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6,2,2),TotVectorTP(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),TotSubjects(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),'VariableNames',{'EMG','TP','Subjects'});
fitlme456Del1 = fitlme(table456Del1,'EMG~TP+(1|Subjects)');
table456Del2 = table(MatrixBinnedREMG(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6,3,2),TotVectorTP(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),TotSubjects(TotVectorTP==4 | TotVectorTP==5 | TotVectorTP==6),'VariableNames',{'EMG','TP','Subjects'});
fitlme456Del2 = fitlme(table456Del2,'EMG~TP+(1|Subjects)');
table789pec1 = table(MatrixBinnedREMG(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9,2,1),TotVectorTP(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9),TotSubjects(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9),'VariableNames',{'EMG','TP','Subjects'});
fitlme789pec1 = fitlme(table789pec1,'EMG~TP+(1|Subjects)');
table789pec2 = table(MatrixBinnedREMG(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9,3,1),TotVectorTP(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9),TotSubjects(TotVectorTP==7 | TotVectorTP==8 | TotVectorTP==9),'VariableNames',{'EMG','TP','Subjects'});
fitlme789pec2 = fitlme(table789pec2,'EMG~TP+(1|Subjects)');
    

%% 3. Figures handles
close all;
% a. Kinematics - general

figure('Name','Kinematics experiment 1','units','normalized','outerposition',[0 0 1 1]); hold on;
subplot(4,10,[1 2 11 12 21 22 31 32]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');% set(gca,'TickLength',[0.025 0.01]);
plot(mean(TotKine(TotVectorTP==1, find(TimeVector==-200):find(TimeVector==500),1),1), mean(TotKine(TotVectorTP==1, find(TimeVector==-200):find(TimeVector==500),2),1),'r','LineWidth',3);
plot(mean(TotKine(TotVectorTP==2, find(TimeVector==-200):find(TimeVector==500),1),1), mean(TotKine(TotVectorTP==2, find(TimeVector==-200):find(TimeVector==500),2),1),'Color',[0 0.4 0],'LineWidth',3);
plot(mean(TotKine(TotVectorTP==3, find(TimeVector==-200):find(TimeVector==500),1),1), mean(TotKine(TotVectorTP==3, find(TimeVector==-200):find(TimeVector==500),2),1),'b','LineWidth',3);
plot(mean(TotKine(TotVectorTP==4, find(TimeVector==-200):find(TimeVector==500),1),1), mean(TotKine(TotVectorTP==4, find(TimeVector==-200):find(TimeVector==500),2),1),'r','LineWidth',3);
plot(mean(TotKine(TotVectorTP==5, find(TimeVector==-200):find(TimeVector==500),1),1), mean(TotKine(TotVectorTP==5, find(TimeVector==-200):find(TimeVector==500),2),1),'Color',[0 0.4 0],'LineWidth',3);
plot(mean(TotKine(TotVectorTP==6, find(TimeVector==-200):find(TimeVector==500),1),1), mean(TotKine(TotVectorTP==6, find(TimeVector==-200):find(TimeVector==500),2),1),'b','LineWidth',3);
plot(mean(TotKine(TotVectorTP==7, find(TimeVector==-200):find(TimeVector==500),1),1), mean(TotKine(TotVectorTP==7, find(TimeVector==-200):find(TimeVector==500),2),1),'r','LineWidth',3);
plot(mean(TotKine(TotVectorTP==8, find(TimeVector==-200):find(TimeVector==500),1),1), mean(TotKine(TotVectorTP==8, find(TimeVector==-200):find(TimeVector==500),2),1),'Color',[0 0.4 0],'LineWidth',3);
plot(mean(TotKine(TotVectorTP==9, find(TimeVector==-200):find(TimeVector==500),1),1), mean(TotKine(TotVectorTP==9, find(TimeVector==-200):find(TimeVector==500),2),1),'b','LineWidth',3);
xlim([-0.06 0.06]); ylim([0.05 0.34]); set(gca,'XColor','none'); set(gca,'YColor','none');
plot(0,0.08,'.','MarkerSize',50); plot(0,0.32,'.','MarkerSize',50);
plot([-0.03 0.03],[0.16 0.16],'k:','LineWidth',3);

subplot(4,10,[7 8 9 10 17 18 19 20]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YAxisLocation','right');
plot(-200:500, mean(TotKineDeviation(TotVectorTPDeviation==4,find(TimeVector==-200):find(TimeVector==500),1),1),'r','LineWidth',3);
plot(-200:500, mean(TotKineDeviation(TotVectorTPDeviation==5,find(TimeVector==-200):find(TimeVector==500),1),1),'Color',[0 0.4 0],'LineWidth',3);
plot(-200:500, mean(TotKineDeviation(TotVectorTPDeviation==6,find(TimeVector==-200):find(TimeVector==500),1),1),'b','LineWidth',3);
plot(-200:500, mean(TotKineDeviation(TotVectorTPDeviation==7,find(TimeVector==-200):find(TimeVector==500),1),1),'r','LineWidth',3);
plot(-200:500, mean(TotKineDeviation(TotVectorTPDeviation==8,find(TimeVector==-200):find(TimeVector==500),1),1),'Color',[0 0.4 0],'LineWidth',3);
plot(-200:500, mean(TotKineDeviation(TotVectorTPDeviation==9,find(TimeVector==-200):find(TimeVector==500),1),1),'b','LineWidth',3);
ylabel('x-deviation [cm]'); yticks([-0.05 -0.025 0 0.025 0.05]); yticklabels({'-5','-2.5','0','2.5','5'}); ylim([-0.06 0.06]);
xlim([-220 520]); xticks([-150 0 150 300 450]); xticklabels({'-0.15','0','0.15','0.3','0.45'}); xlabel('Time [s]');
plot([0 0],[-0.05 0.05],'k:','LineWidth',3);

subplot(4,10,[27 28 29 30]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none'); set(gca,'YAxisLocation','right');
for ii = 1 : nSubjects
    plot([1 2 3], [mean(MatrixHandTrial456(ii,1)), mean(MatrixHandTrial456(ii,2)), mean(MatrixHandTrial456(ii,3))],'.-','Color',[0.7 0.7 0.7],'LineWidth',1.5,'MarkerSize',25);
end
plot([1 2 3], [mean(MatrixHandTrial456(:,1)), mean(MatrixHandTrial456(:,2)), mean(MatrixHandTrial456(:,3))],'k.-','LineWidth',3,'MarkerSize',30);
plot([1 2 3], [0 0 0],'k:','LineWidth',3);
ylim([-0.0075 0.0075]); xlim([-0.5 4.5]);
ylabel('\Delta x-deviation [cm]'); yticks([-0.005 0 0.005]); yticklabels({'-0.5','0','0.5'});

subplot(4,10,[37 38 39 40]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'YAxisLocation','right');
for ii = 1 : nSubjects
    plot([1 2 3], [mean(MatrixHandTrial789(ii,1)), mean(MatrixHandTrial789(ii,2)), mean(MatrixHandTrial789(ii,3))],'.-','Color',[0.7 0.7 0.7],'LineWidth',1.5,'MarkerSize',25);
end
plot([1 2 3], [mean(MatrixHandTrial789(:,1)), mean(MatrixHandTrial789(:,2)), mean(MatrixHandTrial789(:,3))],'k.-','LineWidth',3,'MarkerSize',30);
plot([1 2 3], [0 0 0],'k:','LineWidth',3);
ylim([-0.0075 0.0075]); xlim([-0.5 4.5]);yticks([-0.005 0 0.005]); yticklabels({'-0.5','0','0.5'});
xticks([1 2 3]); xticklabels({'Low','Medium','High'});

subplot(4,10,[3 4 5 6]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');
plot(-200:300,(mean(TotSpeed(TotVectorTP==2,find(TimeVector==-200):find(TimeVector==300),2),1) - mean(TotSpeed(TotVectorTP==1, find(TimeVector==-200):find(TimeVector==300),2),1))./max(mean(TotSpeed(TotVectorTP==1, find(TimeVector==-200):find(TimeVector==300),2),1)),'k','LineWidth',3);
plot(-200:300,(mean(TotSpeed(TotVectorTP==3,find(TimeVector==-200):find(TimeVector==300),2),1) - mean(TotSpeed(TotVectorTP==1, find(TimeVector==-200):find(TimeVector==300),2),1))./max(mean(TotSpeed(TotVectorTP==1, find(TimeVector==-200):find(TimeVector==300),2),1)),'k-.','LineWidth',3);
xlim([-220 320]); ylim([-0.025 0.05]);
yticks([0 0.05]); yticklabels({'0%','5%'});
% yline(0,'k:','LineWidth',2);


subplot(4,10,[13 14 15 16]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');
plot(-200:300,mean([(mean(TotSpeed(TotVectorTP==5,find(TimeVector==-200):find(TimeVector==300),2),1)-mean(TotSpeed(TotVectorTP==4,find(TimeVector==-200):find(TimeVector==300),2),1))./max(mean(TotSpeed(TotVectorTP==4,find(TimeVector==-200):find(TimeVector==300),2),1));(mean(TotSpeed(TotVectorTP==8,find(TimeVector==-200):find(TimeVector==300),2),1)-mean(TotSpeed(TotVectorTP==7,find(TimeVector==-200):find(TimeVector==300),2),1))./max(mean(TotSpeed(TotVectorTP==7,find(TimeVector==-200):find(TimeVector==300),2),1))],1),'k','LineWidth',3);
plot(-200:300,mean([(mean(TotSpeed(TotVectorTP==6,find(TimeVector==-200):find(TimeVector==300),2),1)-mean(TotSpeed(TotVectorTP==4,find(TimeVector==-200):find(TimeVector==300),2),1))./max(mean(TotSpeed(TotVectorTP==4,find(TimeVector==-200):find(TimeVector==300),2),1));(mean(TotSpeed(TotVectorTP==9,find(TimeVector==-200):find(TimeVector==300),2),1)-mean(TotSpeed(TotVectorTP==7,find(TimeVector==-200):find(TimeVector==300),2),1))./max(mean(TotSpeed(TotVectorTP==7,find(TimeVector==-200):find(TimeVector==300),2),1))],1),'k-.','LineWidth',3);
xlim([-220 320]); ylim([-0.02 0.05]);
yticks([0 0.05]); yticklabels({'0%','5%'});
% yline(0,'k:','LineWidth',2);


subplot(4,10,[23 24 25 26]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');
for ii = 1 : nSubjects
    plot([1 2 3],[(MatrixSpeedMeanTrial123(ii,1)), (MatrixSpeedMeanTrial123(ii,2)), (MatrixSpeedMeanTrial123(ii,3))],'.-','Color',[0.7 0.7 0.7],'LineWidth',1.5,'MarkerSize',25);
end
plot([1 2 3],[mean(MatrixSpeedMeanTrial123(:,1)), mean(MatrixSpeedMeanTrial123(:,2)), mean(MatrixSpeedMeanTrial123(:,3))],'k.-','LineWidth',3,'MarkerSize',30);
plot([1 2 3], [0 0 0],'k:','LineWidth',3);
ylim([-0.06 0.06]); xlim([-0.5 4.5]);
ylabel('\Delta speed [cm/s]'); yticks([-0.06 -0.03 0 0.03 0.06]); yticklabels({'-6','-3','0','3','6'});

subplot(4,10,[33 34 35 36]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');
for ii = 1 : nSubjects
    plot([1 2 3],[mean([MatrixSpeedMeanTrial456(ii,1),MatrixSpeedMeanTrial789(ii,1)]), mean([MatrixSpeedMeanTrial456(ii,2),MatrixSpeedMeanTrial789(ii,2)]), mean([MatrixSpeedMeanTrial456(ii,3),MatrixSpeedMeanTrial789(ii,3)])],'.-','Color',[0.7 0.7 0.7],'LineWidth',1.5,'MarkerSize',25);
end
plot([1 2 3],[mean([mean(MatrixSpeedMeanTrial456(:,1)),mean(MatrixSpeedMeanTrial789(:,1))]), mean([mean(MatrixSpeedMeanTrial456(:,2)),mean(MatrixSpeedMeanTrial789(:,2))]), mean([mean(MatrixSpeedMeanTrial456(:,3)),mean(MatrixSpeedMeanTrial789(:,3))])],'k.-','LineWidth',3,'MarkerSize',30);
plot([1 2 3], [0 0 0],'k:','LineWidth',3);
ylim([-0.06 0.06]); xlim([-0.5 4.5]);
xticks([1 2 3]); xticklabels({'Low','Medium','High'});
ylabel('\Delta speed [cm/s]'); yticks([-0.06 -0.03 0 0.03 0.06]); yticklabels({'-6','-3','0','3','6'});


% b. EMG - general 


figure('Name','EMG - Experiment 1','units','normalized','outerposition',[0 0 1 1]); hold on;
subplot(4,8,[1 2 9 10]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');
plot(-200:500, movmean(mean(mean(TotEMG(TotVectorTP==1,find(TimeVector==-200):find(TimeVector==500),:),1),3),5),'r','LineWidth',3);
plot(-200:500, movmean(mean(mean(TotEMG(TotVectorTP==2,find(TimeVector==-200):find(TimeVector==500),:),1),3),5),'Color',[0 0.4 0],'LineWidth',3);
plot(-200:500, movmean(mean(mean(TotEMG(TotVectorTP==3,find(TimeVector==-200):find(TimeVector==500),:),1),3),5),'b','LineWidth',3);
xlim([-120 520]); ylim([0.5 1.2]); yLim = ylim; ylabel('EMG activity [a.u.]');
xline(0,'k:','LineWidth',3)
xline(200,'k:','LineWidth',3);



meandiff = mean(mean(TotEMG(TotVectorTP==1,find(TimeVector==-200):find(TimeVector==500),:),1),3);
subplot(4,8,[17 18 25 26]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');
plot(-200:500, movmean(mean(mean(TotEMG(TotVectorTP==2,find(TimeVector==-200):find(TimeVector==500),:),1),3) - meandiff,10),'k-','LineWidth',3);
plot(-200:500, movmean(mean(mean(TotEMG(TotVectorTP==3,find(TimeVector==-200):find(TimeVector==500),:),1),3) - meandiff,10),'k-.','LineWidth',3);
ylim([-0.1 0.2]); yLim = ylim; xlim([-120 520]); ylabel('EMG activity [a.u.]'); %0-200 vline
xlabel('Time [s]'); xticks([-150 0 200 450]); xticklabels({'-0.15','0','0.2','0.45'});
xline(0,'k:','LineWidth',3); yline(0,'k:','LineWidth',3);
xline(200,'k:','LineWidth',3);

subplot(4,8,[3 4 11 12]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');
for ii = 1 : nSubjects
   plot([1 2 3],MatrixEMGTrial123Pec(ii,:),'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2); 
end
plot([1 2 3],[mean(MatrixEMGTrial123Pec(:,1)), mean(MatrixEMGTrial123Pec(:,2)), mean(MatrixEMGTrial123Pec(:,3))],'k.-','LineWidth',3,'MarkerSize',30);
plot([1 2 3], [0 0 0],'k:','LineWidth',3);
xlim([0 4]); xticks([1 2 3]); xticklabels({'Low','Medium','High'}); ylim([-0.15 0.15]);

% Post hoc analysis

[hph1, pph1,~,sph1] = ttest(MatrixEMGTrial123Pecph(:,1),MatrixEMGTrial123Pecph(:,3),'Tail','left');
% (mean(MatrixEMGTrial123Pecph(:,1))-mean(MatrixEMGTrial123Pecph(:,3)))/std([MatrixEMGTrial123Pecph(:,1);MatrixEMGTrial123Pecph(:,3)])

subplot(4,8,[19 20 27 28]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');
for ii = 1 : nSubjects
    plot([1 2 3],MatrixEMGTrial123Del(ii,:),'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
end
plot([1 2 3],[mean(MatrixEMGTrial123Del(:,1)), mean(MatrixEMGTrial123Del(:,2)), mean(MatrixEMGTrial123Del(:,3))],'k.-','LineWidth',3,'MarkerSize',30);
plot([1 2 3], [0 0 0],'k:','LineWidth',3);
xlim([0 4]); xticks([1 2 3]); xticklabels({'Low','Medium','High'});ylim([-0.25 0.25]);

% Post hoc analysis

[hph2,pph2,~,sph2] = ttest(MatrixEMGTrial123Delph(:,1),MatrixEMGTrial123Delph(:,3),'Tail','left');
% (mean(MatrixEMGTrial123Delph(:,1))-mean(MatrixEMGTrial123Delph(:,3)))/std([MatrixEMGTrial123Delph(:,1);MatrixEMGTrial123Delph(:,3)])

subplot(4,8,[5 6 13 14]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');
plot(-200:500,movmean(mean(mean(TotEMG(TotVectorTP==4,find(TimeVector==-200):find(TimeVector==500),:),1),3),10),'r','LineWidth',3);
plot(-200:500,movmean(mean(mean(TotEMG(TotVectorTP==5,find(TimeVector==-200):find(TimeVector==500),:),1),3),10),'Color',[0 0.4 0],'LineWidth',3);
plot(-200:500,movmean(mean(mean(TotEMG(TotVectorTP==6,find(TimeVector==-200):find(TimeVector==500),:),1),3),10),'b','LineWidth',3);
xlim([-70 220]); xticks([-150 0 50 100 180 450]); xlabel('Time [s]'); xticklabels({'-0.15','0',[],[],[],'0.45'});
ylim([0.5 2]); yLim = ylim;
% plot(linspace(0,0),linspace(yLim(1),yLim(end)),'k:','LineWidth',3);
xline(50,'k:','LineWidth',3); 
xline(100,'k:','LineWidth',3);
xline(180,'k:','LineWidth',3);


meandiff2 = mean(mean(TotEMG(TotVectorTP==4 | TotVectorTP==7,find(TimeVector==-200):find(TimeVector==500),:),1),3);

subplot(4,8,[21 22 29 30]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');
plot(-200:500,movmean(mean(mean(TotEMG(TotVectorTP==8 | TotVectorTP==5,find(TimeVector==-200):find(TimeVector==500),:),1),3) - meandiff2,10),'k-','LineWidth',3);
plot(-200:500,movmean(mean(mean(TotEMG(TotVectorTP==9 | TotVectorTP==6,find(TimeVector==-200):find(TimeVector==500),:),1),3) - meandiff2,10),'k-.','LineWidth',3);
ylim([-0.2 0.5]); yLim = ylim;
xlim([-70 220]); xticks([-150 0 50 100 180 450]); xlabel('Time [s]'); xticklabels({'-0.15','0',[],[],[],'0.45'});
xline(50,'k:','LineWidth',3); yline(0,'k:','LineWidth',3);
xline(100,'k:','LineWidth',3);
xline(180,'k:','LineWidth',3); 

subplot(4,8,[7 8]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none'); set(gca,'YAxisLocation','right');
for ii = 1 : nSubjects
    plot([1 2 3],reshape(MatrixEMGPec456(ii,1,:),1,3),'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
end
plot([1 2 3],[mean(MatrixEMGPec456(:,1,1),1), mean(MatrixEMGPec456(:,1,2),1) mean(MatrixEMGPec456(:,1,3),1)],'k.-','LineWidth',3,'MarkerSize',30); 
plot([1 2 3], [0 0 0],'k:','LineWidth',3);
ylim([-0.3 0.3]); xlim([0 4]);

% Post hoc analysis
[hph3,pph3,~,sph3] = ttest(MatrixEMGPec456ph(:,1,1), MatrixEMGPec456ph(:,1,3),'Tail','left');
% (mean(MatrixEMGPec456ph(:,1,1))-mean(MatrixEMGPec456ph(:,1,3)))/std([MatrixEMGPec456ph(:,1,1);MatrixEMGPec456ph(:,1,3)])

subplot(4,8,[15 16]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');set(gca,'YAxisLocation','right');
for ii = 1 : nSubjects
    plot([1 2 3],reshape(MatrixEMGPec456(ii,2,:),1,3),'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
end
plot([1 2 3],[mean(MatrixEMGPec456(:,2,1),1), mean(MatrixEMGPec456(:,2,2),1) mean(MatrixEMGPec456(:,2,3),1)],'k.-','LineWidth',3,'MarkerSize',30); 
plot([1 2 3], [0 0 0],'k:','LineWidth',3);
ylim([-0.3 0.3]);xlim([0 4]);


%Post hoc analysis 

[hph4,pph4,~,sph4] = ttest(MatrixEMGPec456ph(:,2,1), MatrixEMGPec456ph(:,2,3),'Tail','left');
% (mean(MatrixEMGPec456ph(:,2,1))-mean(MatrixEMGPec456ph(:,2,3)))/std([MatrixEMGPec456ph(:,2,1);MatrixEMGPec456ph(:,2,3)])


subplot(4,8,[23 24]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none'); set(gca,'XColor','none');set(gca,'YAxisLocation','right');
for ii = 1 : nSubjects
    plot([1 2 3],reshape(MatrixEMGDel789(ii,1,:),1,3),'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
end
plot([1 2 3], [mean(MatrixEMGDel789(:,1,1),1), mean(MatrixEMGDel789(:,1,2),1), mean(MatrixEMGDel789(:,1,3),1)],'k.-','LineWidth',3,'MarkerSize',30);
plot([1 2 3], [0 0 0],'k:','LineWidth',3);
ylim([-0.6 0.6]);xlim([0 4]);

% Post hoc analysis 

[hph5,pph5,~,sph5] = ttest(MatrixEMGDel789ph(:,1,1), MatrixEMGDel789ph(:,1,3),'Tail','left');
% (mean(MatrixEMGDel789ph(:,1,1))-mean(MatrixEMGDel789ph(:,1,3)))/std([MatrixEMGDel789ph(:,1,1);MatrixEMGDel789ph(:,1,3)])

subplot(4,8,[31 32]); hold on; set(gca,'LineWidth',4); set(gca,'FontSize',18); set(gca,'Color','none');set(gca,'YAxisLocation','right');
for ii = 1 : nSubjects
    plot([1 2 3],reshape(MatrixEMGDel789(ii,2,:),1,3),'.-','Color',[0.7 0.7 0.7],'MarkerSize',25,'LineWidth',2);
end
plot([1 2 3], [mean(MatrixEMGDel789(:,2,1),1), mean(MatrixEMGDel789(:,2,2),1), mean(MatrixEMGDel789(:,2,3),1)],'k.-','LineWidth',3,'MarkerSize',30);
plot([1 2 3], [0 0 0],'k:','LineWidth',3);
xticks([1 2 3]); xticklabels({'Low','Medium','High'}); ylim([-0.6 0.6]);xlim([0 4]);

% Post hoc analysis 

[hph6,pph6,~,sph6] = ttest(MatrixEMGDel789ph(:,2,1), MatrixEMGDel789ph(:,2,3),'Tail','left');
% (mean(MatrixEMGDel789ph(:,2,1))-mean(MatrixEMGDel789ph(:,2,3)))/std([MatrixEMGDel789ph(:,2,1);MatrixEMGDel789ph(:,2,3)])
