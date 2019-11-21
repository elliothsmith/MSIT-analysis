function [coherenceStats] = analyzeMSITcoherence(patientID,sessionNum,nevFile)
%ANALYZEMSITCOHERENCE does spike-field coherence analysis on MSIT data.
%
%   [coherenceStats] = analyzeMSITcoherence(patientID,sessionNum,nevFile)
%   calculates spike-field coherograms for MSIT data in nevFile and its
%   associated ns(3) file. analyzeMSITcoherence saves statistics, and
%   plots results.
%
%   Two tips for easy usage:
%       1) place the nev file in the same directory as the nsx files.
%       2) make sure that the nev file with sorted units has its original
%           name at the beginning of the file name.
%


% author: ElliotHSmith (https://github.com/elliothsmith/MSIT-analysis)



%% loading data from NEV file
display('loading action potential and local field potential data...')
[dataPath, nvName, nvExt] = fileparts(nevFile);

% defining nsFile.
try
    nsFile = fullfile(dataPath,[nvName(1:19) '.ns3']);
catch
    nsFile = fullfile(dataPath,[nvName(1:8) '*.ns3']);
end


% parsing files and loading data.
if strcmp(nvExt,'.nev')
    NEV = openNEV(nevFile,'read');
    if ~exist(nsFile,'file')
        [nsFile,nsPath,~] = uigetfile('*.ns*','Select the correct NSx file');
        NS3 = openNSx(fullfile(nsPath,nsFile));
    else
        NS3 = openNSx(nsFile);
    end
elseif strcmp(nvExt,'.mat')
    load(nevFile);
    if ~exist(nsFile,'file')
        [nsFile,nsPath,~] = uigetfile('*.ns*','Select the correct NSx file');
        NS3 = openNSx(fullfile(nsPath,nsFile));
    else
        NS3 = openNSx(nsFile);
    end
end


%% organizing important task parameters.
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% [20160824] removing practice trial triggers and updating number of trials.
if strcmp(patientID,'CUBF09')
    trigs(95) = 104;
    [~,trigTimes] = removePracticeTriggers(trigs,trigTimes);
else
    [trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
end
nTrials = sum(trigs==90);


%% discharge detection GUI.
rejectTrials = dischargeRejectionTool(patientID,sessionNum,nevFile,0);


%% organizing responses.
responses = trigs(trigs>=100 & trigs<=104);
cues = trigs(trigs>=1 & trigs<28);
noResponseIdcs = responses==100;


%% reaction time calculation
rt = responses-cues;


%% defining Chronux parameters.
movingWin = [1 0.010];
params.Fs = 2e3; % sampling frequency for LFP
params.pad = 1; %
params.fpass = [0 50]; % frequency range of interest
params.tapers = [5 9]; % emphasize smoothing for the spikes
params.trialave = 1; % average over trials {CHANGES BELOW}
params.err = [1 0.05]; % population error bars


%% creating neural timing variable
if strcmp(patientID,'CUBF09')
    ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode) double(NEV.Data.Spikes.Unit) double(NEV.Data.Spikes.TimeStamp)];
else
    ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];
end
ChanUnitTimestamp(ChanUnitTimestamp(:,2)==255,:) = [];


%% parsing channel labels
chanswUnits = unique(ChanUnitTimestamp(:,1));
nChans = length(chanswUnits);
unitChanLabels = [NEV.ElectrodesInfo.ElectrodeLabel]';
inclMicroChannels = unitChanLabels(unitChanLabels(:,1)==repmat('u',size(unitChanLabels(:,1))),:);
numBFs = floor(size(inclMicroChannels,1)./8); % hard code justification: 8 is the number of BF microcontacts
for bf = 1:numBFs
    tmp = deblank(inclMicroChannels((bf)*8,:));
    BFlabels{bf} = tmp(1:end-1);
end
BFlabels


%% aligning data on specific or all trial markers: query user.
% alignSpot = input('\nAlign on: \n1)fixation? \n2)cue? \n3)response? \n4)all of the above ?');
% if isequal(alignSpot,4)
%     aS = 1:3;
% else
%     aS = alignSpot;
% end


display('Aligning AP data on stimulus and response.');
% for loop to save multiple epochs
for aS = 1
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            nTrials = length(trialStarts);
            pre = 2;
            post = 5;
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>100 & trigs<=105);
            nTrials = length(trialStarts);
            pre = 3;
            post = 5;
    end    
    % looping over Channels
    for ch = 1:nChans
        % including labels in the data structure - doesn't mess up chronux
        data(aS).channel(ch).label = deblank(unitChanLabels(chanswUnits(ch),:));
        % looping over number of units in the AP data
        nUnits = length(unique(ChanUnitTimestamp(chanswUnits(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
        for un = 1:nUnits
            % getting unit times for the current channel and unit.
            unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==chanswUnits(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds
            % loooping over trials
            for tt = 1:nTrials
                %% putting the data in a structure.
                data(aS).channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>(trialStarts(tt)-pre) & unitTimes<trialStarts(tt)+post)...
                    - repmat(trialStarts(tt)-pre,length(sum(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post)),1);
            end
        end
    end
    % save data structure?
    saveFlag = 0;
    if saveFlag
        display(['saving spike time structure for ' alignName '...'])
        save([patientID '_session' num2str(sessionNum) '_spikeTimeStruct_alignedon' alignName '.mat'],'data','pre','post','alignName')
    end
end


%%  alinging LFP data on the same time base as AP data.
tmp = double(NS3.Data);
LFPlabels = {NS3.ElectrodesInfo.Label};
%% TODO: figure out a smart way of allocating which LFP contacts you want
display('Aligning LFP data on stimulus and response.');
% for loop to save multiple epochs
for aS = 1:2
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            nTrials = length(trialStarts);
            pre = 2;
            post = 5;
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>100 & trigs<=105);
            nTrials = length(trialStarts);
            pre = 3;
            post = 5;
    end
    
    % initializing LFPmat
    LFPmat{aS} = zeros(((pre+post)*params.Fs)+1,nTrials,length(LFPlabels));
    
    
    %% looping over channels and trials
    for ch = 1:length(LFPlabels)
        for tt = 1:nTrials
            %% tensorizing:
            try
                LFPmat{aS}(:,tt,ch) = tmp(ch,floor((trialStarts(tt)-pre)*params.Fs):floor((trialStarts(tt)+post)*params.Fs));
            catch
                LFPmat{aS}(:,tt,ch) = tmp(ch,floor((trialStarts(tt)-pre)*params.Fs):floor((trialStarts(tt)-pre)*params.Fs)+(size(LFPmat{aS},1)-1));
            end
            
        end
    end
end


%% parsing behavior
trialType = zeros(1,nTrials);
condition = trigs(trigs>=1 & trigs<=27);


%% setting up codes for PSTHs over conflict types.
% These are the correct codes. Double Checked on 20160216
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)


%% Rasters and PSTHs
for aS2 = 1:length(data)
    % which alignment spot
    switch aS2
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            nTrials = length(trialStarts);
            pre = 2;
            post = 5;
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>100 & trigs<=105);
            nTrials = length(trialStarts);
            pre = 3;
            post = 5;
    end
    
    
    for ch = 1:size(data(aS2).channel,2)
        % [20160607] making sure that spike field coherence is calculated
        % between a spike and its local and distal LFPs on the same BF
        tmp = deblank(data(aS2).channel(ch).label);
        for bf = 1:length(BFlabels)
            if strcmp(tmp(1:end-1),BFlabels{bf})
                contChans = [1 8] + repmat((bf-1)*8,1,2); % hard code justification: 8 is the number of BF microcontacts
            end
        end
        for un = 1:size(data(aS2).channel(ch).unit,2)
            for cCH = contChans
                % [20160607] setting up labels
                if mod(cCH,8)
                    LFPcohLoc = 'Medial';
                else
                    LFPcohLoc = 'Lateral';
                end
                
                
                %% [20160622] rejecting the set comprising the union of the bad trials from all channels.
                rejected = union(rejectTrials(contChans(1)).rejectTheseTrials,rejectTrials(contChans(2)).rejectTheseTrials);
                notRejected = ones(1,nTrials);
                notRejected(rejected) = false;
                trialType = trialType.*notRejected;
                trialTypeRJCTD = trialType(logical(notRejected));
                
                
                if isequal(params.trialave,0)
                    %% calculating coherogram for each trial
                    %                 display(sprintf('calculating spike-field coherence within trials for\n    LFP channel %s,\n    AP channel %s,\n    unit %d,\n    aligned on %s.'...
                    %                     ,strtrim(LFPlabels{cCH}),deblank(unitChanLabels(ch,1:5)),un,alignName));
                    
                    % caluclate coherence without averaging over trials.
                    tic
                    [C_alltrials,phi_alltrials,~,~,~,t,f] = cohgramcpt(squeeze(LFPmat{aS2}(:,:,cCH)),data(aS2).channel(ch).unit(un).trial, movingWin, params);
                    tspec = linspace(t(1)-pre,t(end)-pre,length(t));
                    toc
                    
%                     [nullDist,nullPhase] = generateCoherenceNullDist(500, squeeze(LFPmat{aS2}(:,:,cCH)), data(aS2).channel(ch).unit(un).trial, movingWin, params);

                    

                elseif isequal(params.trialave,1)
                    %% calculating avraged cohereograms for each condition
                    %                     display(sprintf('calculating spike-field coherence across trials for\n    LFP channel %s,\n    AP channel %s,\n    unit %d,\n    aligned on %s.'...
                    %                         ,LFPlabels,strtrim(unitChanLabels(ch,1:5)),un,alignName));
                    % coherence calculation for each trialType
                    
                    savedFileName = sprintf('/media/user1/data4TB/msit_units/Experiment_I__ACC/%s/Data/coherence/%sAlignedCoherence_averagedCohgrams_session%d_unit%d_channel%d_continuousChannel%s.mat',patientID,alignName,sessionNum,un,ch,LFPcohLoc);
                    %
                    %                 if exist(savedFileName,'file')
                    %                     load(savedFileName)
                    %                 else
                    
                    
                    timeYouWantThisToTake = 'way long'; % flag for calculating permutation distributions
                    if strcmp(timeYouWantThisToTake,'way long')
%                         tic
%                         display('calculating trial-averaged coherence for easy trials...')
%                         [C0,phi0,~,~,~,t,f,~,phiStd0,Cerr0] ...
%                             = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==1,cCH)),data(aS2).channel(ch).unit(un).trial(trialType==1), movingWin, params);
%                         display('calculating trial-averaged coherence for simon trials...')
%                         [C1a,phi1a,~,~,~,t,f,~,phiStd1a,Cerr1a] ...
%                             = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==2,cCH)),data(aS2).channel(ch).unit(un).trial(trialType==2), movingWin, params);
%                         display('calculating trial-averaged coherence for flanker trials...')
%                         [C1b,phi1b,~,~,~,t,f,~,phiStd1b,Cerr1b] ...
%                             = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==3,cCH)),data(aS2).channel(ch).unit(un).trial(trialType==3), movingWin, params);
%                         display('calculating trial-averaged coherence for hard trials...')
%                         [C2,phi2,~,~,~,t,f,~,phiStd2,Cerr2] ...
%                             = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==4,cCH)),data(aS2).channel(ch).unit(un).trial(trialType==4), movingWin, params);
%                          tspec = linspace(t(1)-pre,t(end)-pre,length(t));
%                         toc
                        
                        
                        %% [20180529]::revisions:: controlling for numbers of trials. 
                        cond0 = find(trialType==1);
                        cond1a = find(trialType==2);
                        cond1b = find(trialType==3);
                        cond2 = find(trialType==4);
                        nTrials_controlled = min([length(cond0) length(cond1a) length(cond1b) length(cond2)]);
                        cond0 = datasample(cond0,nTrials_controlled);
                        cond1a = datasample(cond1a,nTrials_controlled);
                        cond1b = datasample(cond1b,nTrials_controlled);
                        cond2 = datasample(cond2,nTrials_controlled);
                        
                        
                        tic
                        display('calculating trial-averaged coherence for easy trials...')
                        [C0,phi0,~,~,~,t,f,~,phiStd0,Cerr0] ...
                            = cohgramcpt(squeeze(LFPmat{aS2}(:,cond0,cCH)),data(aS2).channel(ch).unit(un).trial(cond0), movingWin, params);
                        display('calculating trial-averaged coherence for simon trials...')
                        [C1a,phi1a,~,~,~,t,f,~,phiStd1a,Cerr1a] ...
                            = cohgramcpt(squeeze(LFPmat{aS2}(:,cond1a,cCH)),data(aS2).channel(ch).unit(un).trial(cond1a), movingWin, params);
                        display('calculating trial-averaged coherence for flanker trials...')
                        [C1b,phi1b,~,~,~,t,f,~,phiStd1b,Cerr1b] ...
                            = cohgramcpt(squeeze(LFPmat{aS2}(:,cond1b,cCH)),data(aS2).channel(ch).unit(un).trial(cond1b), movingWin, params);
                        display('calculating trial-averaged coherence for hard trials...')
                        [C2,phi2,~,~,~,t,f,~,phiStd2,Cerr2] ...
                            = cohgramcpt(squeeze(LFPmat{aS2}(:,cond2,cCH)),data(aS2).channel(ch).unit(un).trial(cond2), movingWin, params);
                        tspec = linspace(t(1)-pre,t(end)-pre,length(t));
                        toc
                        
                        
                        %% [20160920] saving coherence data over trials.
                        saveData=1;
                        if saveData
                            save(savedFileName,'C0','C1a','C1b','C2',...
                                'phi0','phi1a','phi1b','phi2',...
                                'phiStd0','phiStd1a','phiStd1b','phiStd2',...
                                'Cerr0','Cerr1a','Cerr1b','Cerr2','t','f',...
                                'params','movingWin','patientID','tspec','-v7.3')    %,'nullDist','nullPhase'
                        end
                        
                        statAlpha = mean([min(Cerr0) min(Cerr1a) min(Cerr1b) min(Cerr2)]);
                        
                    elseif strcmp(timeYouWantThisToTake,'short please')
                        tic
                        display('calculating trial-averaged coherence for easy trials...')
                        [C0,phi0,~,~,~,~,f] ...
                            = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==1,cCH)),data(aS2).channel(ch).unit(un).trial(trialType==1), movingWin, params);
                        display('calculating trial-averaged coherence for simon trials...')
                        [C1a,phi1a,~,~,~,t,f] ...
                            = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==2,cCH)),data(aS2).channel(ch).unit(un).trial(trialType==2), movingWin, params);
                        display('calculating trial-averaged coherence for flanker trials...')
                        [C1b,phi1b,~,~,~,t,f] ...
                            = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==3,cCH)),data(aS2).channel(ch).unit(un).trial(trialType==3), movingWin, params);
                        display('calculating trial-averaged coherence for hard trials...')
                        [C2,phi2,~,~,~,t,f] ...
                            = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==4,cCH)),data(aS2).channel(ch).unit(un).trial(trialType==4), movingWin, params);
                        
                        tspec = linspace(t(1)-pre,t(end)-pre,length(t));
                        toc
                        
                        % just for plotting!!
                        statAlpha = 0.05;
                        Cerr0  = zeros(size(C0));
                        Cerr1a  = zeros(size(C1a));
                        Cerr1b  = zeros(size(C1b));
                        Cerr2  = zeros(size(C2));
                    end

                end
                
%                 % saving
% %                 maximize(ch*1000+un)
%                 fName = sprintf('./Figs/coherence/%s_session_%d_Channel_%d_Unit_%d_Coherence_with%sLFP_%saligned',patientID,sessionNum,chanswUnits(ch),un,LFPcohLoc,alignName);
%                 print(ch*1000+un,fName, '-dpdf','-bestfit')
%                 close(ch*1000+un)
                

            end % looping over lfp channels
        end % looping over units for each channel and align spot
    end % looping over channels for each align spot.
    clear trialtype notRejected
end % looping over align spots (Stimulus & response)

