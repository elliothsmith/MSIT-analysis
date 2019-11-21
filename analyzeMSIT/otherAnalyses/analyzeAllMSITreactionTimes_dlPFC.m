% parse behavior take a BHV mat file from read_bhv and parse the trial
% parameters [20150109]
clear; clc; close all;


%% //conflict colors//
col0 = [183 30 103]./255;
col1a = [246 139 31]./255;
col1b = [0 166 81]./255;
col2 = [82 79 161]./255;


RTchoice = 1 %input('raw RTs (1) , z-scored RTs (2) , or speed of target selection (3) ?');


%% [20170303] initializing variables. 
allRTs =[];
allConditions = [];
subjectSessions = [];
subjects = [];

for pt = 1:8
    switch pt
        case 1
            patientNumber = 10;
%         case 2
%             patientNumber = 17;
        case 2
            patientNumber = 18;
        case 3
            patientNumber = 20;
        case 4
            patientNumber = 21;
        case 5
            patientNumber = 23;
        case 6
            patientNumber = 26;
        case 7
            patientNumber = 30;
        case 8
            patientNumber = 33;
    end
    
    display('loading data...')
    patient = ['C' num2str(patientNumber)];
    
    % % load file from read_bhv
    % load('C:\Users\Elliot\Documents\MATLAB\Analysis\MSIT\dlPFC_OR\data\behav_one.mat')
    % bhv = struct();
    
    cd(['~/data/msit_units/Experiment_II__dlPFC/' patient])
    
    if strcmp(patient,'C10')
        % load file for c17
        load('behav_one.mat')
        bhv = struct();
        startTrial = 1;
    elseif strcmp(patient,'C17')
        % load file for c17
        load('******21_run1.mat')
        behav_one = BHV;
        bhv = struct();
        startTrial = 1;
    elseif strcmp(patient,'C18')
        % load file for c18
        load('******11_run1.mat')
        behav_one = BHV;
        bhv = struct();
        startTrial = 23;
    elseif strcmp(patient,'C20')
        % load file for c20
        load('c20_ORMSIT_BHV.mat')
        behav_one = BHV;
        bhv = struct();
        startTrial = 10; % This is the trial the recording started on
    elseif strcmp(patient,'C21')
        % load file for c21
        load('c21_ORMSIT_BHV.mat')
        behav_one = BHV;
        bhv = struct();
        startTrial = 1;
    elseif strcmp(patient,'C23')
        % load file for c21
        load('c23_bhv.mat')
        behav_one = BHV;
        bhv = struct();
        startTrial = 21;
    elseif strcmp(patient,'C26')
        % load file for c21
        load('c26_bhv.mat')
        behav_one = bhv;
        bhv = struct();
        startTrial = 57;
    elseif strcmp(patient,'C30')
        % load file for c30
        load('C30_MSIT_BHV.mat')
        behav_one = bhv;
        bhv = struct();
        startTrial = 1;
    elseif strcmp(patient,'C33')
        % load file for c30
        load('C33_MSIT_OR_LFP_BHV.mat')
        behav_one = bhv;
        bhv = struct();
        startTrial = 1;
    end
    
    % saving variables of interest
    bhv.file = behav_one.DataFileName;
    bhv.startTime = behav_one.StartTime;
    bhv.endTime = behav_one.FinishTime;
    
    % saving reaction times
    bhv.rt = behav_one.ReactionTime(startTrial:end);
    
    length(bhv.rt)
    
    timeouts = isnan(bhv.rt);
    bhv.rt(timeouts)=[];
    
    % saving error trials.
    bhv.errorTrials = behav_one.TrialError(startTrial:end);
    bhv.errorTrials(timeouts)=[];
    
    % easy medium and hard trials.
    bhv.conditions = behav_one.ConditionNumber(startTrial:end);
    bhv.conditions(timeouts)=[];
    bhv.type0 = logical(bhv.conditions>=1 & bhv.conditions<4);   % Type 0 (Cond # 1-3)
    bhv.type2 = logical(bhv.conditions>3 & bhv.conditions<16);   % Type 2 (Cond # 4-15)
    bhv.type1a = logical(bhv.conditions>15 & bhv.conditions<19);   % Type 1a Spatial interference (Cond # 16-18)
    bhv.type1b = logical(bhv.conditions>18 & bhv.conditions<28);   % Type 1b Distractor interference (Cond # 16-27)
    conditions = bhv.type0+(2*bhv.type1a)+(3*bhv.type1b)+(4*bhv.type2);
    
    
    %% [20170303] Finally doing things right and setting up a mixed effects model
    allRTs = cat(2,allRTs,bhv.rt);
    allConditions = cat(2,allConditions,conditions');
    subjectSessions = cat(2,subjectSessions,repmat(pt,1,length(bhv.rt)));
    
    
    %% [20160319] setting up indices.
    % conflict indices fro previous trials.
    noConfTrials = bhv.type0;
    allConfTrials = bhv.type1a + bhv.type1b + bhv.type2;
    hardTrials = bhv.type2;
    
    
    %% post-error behavior
    postError = logical([bhv.errorTrials(2:end); 0]);
    statusQuo = [bhv.errorTrials(2:end); 0]==0;
    
    
    %% feedback behavior
    %     postCorrectFB
    %     postErrorFB
    %     postNeutralFB
    %
    %     postTransValenced
    %     postTransNeutral
    
    
    %% which type of RT do you like?
    switch RTchoice
        case 1
            RTtype = 'rawRTs';
            
            %% behavior by conflict
            RTs0{pt} = bhv.rt(bhv.type0);
            RTs1a{pt} = bhv.rt(bhv.type1a);
            RTs1b{pt} = bhv.rt(bhv.type1b);
            RTs2{pt} = bhv.rt(bhv.type2);
            
%             % visualizing each RT distribution
%             figure(pt)
%             hold on
%             h0 = histfit(RTs0{pt});
%             h0(1).EdgeColor = col0;
%             h0(1).FaceColor = 'none';
%             h0(2).Color = col0;
%             h1a = histfit(RTs1a{pt});
%             h1a(1).EdgeColor = col1a;
%             h1a(1).FaceColor = 'none';
%             h1a(2).Color = col1a;
%             h1b = histfit(RTs1b{pt});
%             h1b(1).EdgeColor = col1b;
%             h1b(1).FaceColor = 'none';
%             h1b(2).Color = col1b;
%             h2 = histfit(RTs2{pt});
%             h2(1).EdgeColor = col2;
%             h2(1).FaceColor = 'none';
%             h2(2).Color = col2;
%             hold off
%             pause
            
            
            %% behavior by previous trial effects.
            % previous trial no conflict
            RTs0b40{pt} = bhv.rt([0; diff(allConfTrials)]==0 & noConfTrials==1);
            RTs0b42{pt} = bhv.rt([0; diff(hardTrials)]==-1);
            
            % previous trial conflict
            RTs12b40{pt} = bhv.rt([0; diff(allConfTrials)]==1);
            RTs12b42{pt} = bhv.rt([0; diff(allConfTrials)]==0 & hardTrials);
            
            % removing the first trial, since there isn't a trial before it
            RTs0b40{pt}(1) = [];
            RTs0b42{pt}(1) = [];
            RTs12b40{pt}(1) = [];
            RTs12b42{pt}(1) = [];
            
            %% post error behavior
            postErrorRTs{pt} = bhv.rt(postError);
            statusQuoRTs{pt} = bhv.rt(statusQuo);
            
            %% post-feedback behavior.
            
            
        case 2
            RTtype = 'zScoredRTs';
            
            %% [20160319] z-scored RTs
            % and z-scored RTs
            RTs0{pt} = zscore(bhv.rt(bhv.type0));
            RTs1a{pt} = zscore(bhv.rt(bhv.type1a));
            RTs1b{pt} = zscore(bhv.rt(bhv.type1b));
            RTs2{pt} = zscore(bhv.rt(bhv.type2));
            
            % previous trial no conflict
            RTs0b40{pt} = zscore(bhv.rt([0; diff(allConfTrials)]==0 & noConfTrials==1));
            RTs0b42{pt} = zscore(bhv.rt([0; diff(hardTrials)]==-1));
            
            % previous trial conflict
            RTs12b40{pt} = zscore(bhv.rt([0; diff(allConfTrials)]==1));
            RTs12b42{pt} = zscore(bhv.rt([0; diff(allConfTrials)]==0 & hardTrials));
            
            % removing the first trial, since there isn't a trial before it
            RTs0b40{pt}(1) = [];
            RTs0b42{pt}(1) = [];
            RTs12b40{pt}(1) = [];
            RTs12b42{pt}(1) = [];
            
            %% post-error behavior
            postErrorRTs{pt} = zscore(bhv.rt(postError));
            statusQuoRTs{pt} = zscore(bhv.rt(statusQuo));
            
            %% post-feedback behavior.
            
            
        case 3
            RTtype = 'STS';
            
            %% [20160319] speed of target selection
            RTbar = nanmean(bhv.rt);
            
            % normalized STSs
            RTs0{pt} = (1./bhv.rt(bhv.type0))./RTbar;
            RTs1a{pt} = (1./bhv.rt(bhv.type1a))./RTbar;
            RTs1b{pt} = (1./bhv.rt(bhv.type1b))./RTbar;
            RTs2{pt} = (1./bhv.rt(bhv.type2))./RTbar;
            
            % previous trial no conflict
            RTs0b40{pt} = (1./bhv.rt([0; diff(allConfTrials)]==0 & noConfTrials==1))./RTbar;
            RTs0b42{pt} = (1./bhv.rt([0; diff(hardTrials)]==-1))./RTbar;
            
            % previous trial conflict
            RTs12b40{pt} = (1./bhv.rt([0; diff(allConfTrials)]==1))./RTbar;
            RTs12b42{pt} = (1./bhv.rt([0; diff(allConfTrials)]==0 & hardTrials))./RTbar;
            
            % removing the first trials
            RTs0b40{pt}(1) = [];
            RTs0b42{pt}(1) = [];
            RTs12b40{pt}(1) = [];
            RTs12b42{pt}(1) = [];
            
            %% post-error behavior
            postErrorRTs{pt} = (1./bhv.rt(postError))/RTbar;
            statusQuoRTs{pt} = (1./bhv.rt(statusQuo))/RTbar;
            
            %% post-feedback behavior.
            
            
            
            
    end
end


%% [20161103] now adding data from Utah Array pt.
load('/home/user1/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/CUCX2_behavioralStatistics_OverSessions.mat')
% the last patient. 
pt = length(RTs0)+1;

% % adding cereplex patient to RT distributions.\
subjects = cat(2,subjectSessions,repmat(subjectSessions(end)+1,1,length([RTs{1}])));
allRTs = cat(2,allRTs./1000,RTs{1});
allConditions = cat(2,allConditions,trialType{1});
subjectSessions = cat(2,subjectSessions,repmat(9,1,length(RTs{1})));

% removing jitter timeout trials. 
jitterTimeouts = [allRTs>4.5 & allRTs<0.2];
allConditions(jitterTimeouts) = [];
allConditions(allConditions==1) = 'A';
allConditions(allConditions==2) = 'B';
allConditions(allConditions==3) = 'C';
allConditions(allConditions==4) = 'D';
allRTs(jitterTimeouts) = [];
subjectSessions(jitterTimeouts) = [];
subjects(jitterTimeouts) = [];


% %% [20180503] plot all subject RTs. 
% for fck = 1:length(unique(subjectSessions))
% figure(99)
% hold on
% violinPlot(allRTs(subjectSessions==fck & allConditions=='A')','showMM',4,'color',col0,'xValues',(fck-1)*4)
% violinPlot(allRTs(subjectSessions==fck & allConditions=='B')','showMM',4,'color',col1a,'xValues',((fck-1)*4)+1)
% violinPlot(allRTs(subjectSessions==fck & allConditions=='C')','showMM',4,'color',col1b,'xValues',((fck-1)*4)+2)
% violinPlot(allRTs(subjectSessions==fck & allConditions=='D')','showMM',4,'color',col2,'xValues',((fck-1)*4)+3)
% 
% end
% hold off
% maximize(99)
% ylim([0 5])
% orient(99,'landscape')
% saveas(99,'~/Dropbox/allRTs_allSessions_dlPFC.pdf')

precentLessThan50 = (sum(allRTs>0 & allRTs<0.5)./length(allRTs))*100;
precentLessThan75 = (sum(allRTs>0 & allRTs<0.75)./length(allRTs))*100;


%% [20170303] running the mixed effects model. 
RTtbl = table(log10(allRTs)',allConditions',subjectSessions',subjects','VariableNames',{'RT','condition','session','subject'});
glme = fitglme(RTtbl,'RT ~ 1 + condition*session + (1 + condition*session|subject)')
plotResiduals(glme,'probability','ResidualType','Pearson')
[h,p,k,c] = lillietest(real(log(glme.residuals)))


