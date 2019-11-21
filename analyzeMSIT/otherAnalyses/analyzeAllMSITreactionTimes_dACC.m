% [20160822] analyze all reaction times for experiment I
clear all

%% //conflict colors//
col0 = [183 30 103]./255;
col1a = [246 139 31]./255;
col1b = [0 166 81]./255;
col2 = [82 79 161]./255;
confColors = [col0; col1a; col1b; col2];

% finding beghavioral statistics from analyzeMSITbehavior.m
behavioralDataList = subdir('~/data/msit_units/*behavioralStatistics_OverSessions.mat');

% initializing vector of RTs.
allRTs = []; allTrialTypes = []; allErrors = []; subjectSessions = []; sub = 0; perSessErrorRate = []; subjects = [];
for fl = 1:length(behavioralDataList) % looping over each patient
    if isempty(strfind(behavioralDataList(fl).name,'CUCX2')) % excluding CUCX2
        %         behavioralDataList(fl).name(1:6)
        % loading data from each patient
        load(behavioralDataList(fl).name)
        behavioralDataList(fl).name
        for s = 1:length(RTs) % looping over sessions per patient
            % putting each RT in a verctor of all RTs
            allRTs = cat(2,allRTs,(RTs{s}));
            allTrialTypes = cat(2,allTrialTypes,trialType{s}(1:length(RTs{s})));
            allErrors = cat(1,allErrors,errors{s});
            sessionIdx = sub+s;
            length(RTs{s})
            subjectSessions = cat(2,subjectSessions,repmat(sessionIdx,1,length(RTs{s})));
            %             perSessErrorRate = cat(2,perSessErrorRate,nErrors{s}./length(RTs{s}));
            subjects = cat(2,subjects,repmat(fl,1,length(RTs{s})));
            
        end
        sub = subjectSessions(end);
    end
    
    for cnd = 1:4
        meanRTs(cnd,fl) = nanmedian(allRTs(allTrialTypes==cnd & subjects==fl));
        stdRTs(cnd,fl) = nanstd(allRTs(allTrialTypes==cnd & subjects==fl));
    end
    
end

% removing jitter timeout trials.
jitterTimeouts = (allRTs>=4.5 | allRTs<0.2);
allTrialTypes(jitterTimeouts) = [];
allTrialTypes(allTrialTypes==1) = 'A';
allTrialTypes(allTrialTypes==2) = 'B';
allTrialTypes(allTrialTypes==3) = 'C';
allTrialTypes(allTrialTypes==4) = 'D';
allRTs(jitterTimeouts) = [];
subjectSessions(jitterTimeouts) = [];
subjects(jitterTimeouts) = [];
allErrors(jitterTimeouts) = [];

% LMMs
RTtbl = table((allRTs(subjects<6))',allTrialTypes(subjects<6)',subjectSessions(subjects<6)',subjects(subjects<6)','VariableNames',{'RT','condition','session','subject'});
glme = fitglme(RTtbl,'RT ~ 1 + condition*session + (1 + condition*session | subject)')
plotResiduals(glme,'probability','ResidualType','Pearson')
[h,p,k,c] = lillietest(glme.residuals)


