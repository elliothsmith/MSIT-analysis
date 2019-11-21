function [coherenceStats] = analyzeMSITpopulationCoherenceU(patientID,sessionNum,nevFile)
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
nsFile = fullfile(dataPath,[nvName(1:19) '.ns3']);
if ~exist(nsFile,'file')
    nsList = dir(fullfile(dataPath,[nvName(1:8) '*.ns3']));
    nsFile = fullfile(dataPath,nsList(1).name);
end


% if nsFile == -1;
%     display('could not find ECoG data. Loading full bandwidth data...')
%     nsFile = fullfile(dataPath,[nvName(1:19) '.ns5']);
%     utahArrayFlag = 1;
% end

% parsing files and loading data.
if strcmp(nvExt,'.nev')
    NEV = openNEV(nevFile,'read');
elseif strcmp(nvExt,'.mat')
    load(nevFile);
end

% loading LFP data.
if ~exist(nsFile,'file')
    [nsFile,nsPath,~] = uigetfile('*.ns*','Select the correct NSx file');
    NS3 = openNSx(fullfile(nsPath,nsFile));
else
    NS3 = openNSx(nsFile);
end


%% organizing important task parameters.
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% [20160824] removing practice trial triggers and updating number of trials.
if ~(strcmp(patientID,'CUCX2') && isequal(sessionNum,3))
    [trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
    nTrials = sum(trigs==90);
end


%% [20170329] removing timed out trials
timeOutTrials = removeTimeoutTrials(trigs);


%% discharge detection GUI.
rejectTrials = dischargeRejectionToolU(patientID,sessionNum,nevFile,0);


%% organizing responses.
responses = trigs(trigs>100 & trigs<=104);


%% reaction time calculation
rt = trigTimes(trigs>=100 & trigs<104) - trigTimes(trigs>=1 & trigs<28);


%% defining Chronux parameters.
movingWin = [1 0.025];
params.Fs = 2e3; % sampling frequency for LFP
params.pad = 0; %
params.fpass = [0 50]; % frequency range of interest
params.tapers = [5 9]; % emphasize smoothing for the spikes
params.trialave = 1; % average over trials {CHANGES BELOW}
params.err = [1 0.01]; % population error bars


%% creating neural timing variable
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];
ChanUnitTimestamp(ChanUnitTimestamp(:,2)==255,:) = [];


%% parsing channel labels
chanswUnits = unique(ChanUnitTimestamp(:,1));
nChans = length(chanswUnits);
unitChanLabels = [NEV.ElectrodesInfo.ElectrodeLabel]';


%% calculating the total number of units.
for ch = 1:nChans
    ChanUnits(ch) = max(unique(ChanUnitTimestamp(ChanUnitTimestamp(:,1)==chanswUnits(ch),2)));
    % creating a little array if channel and unit indices for each "Neuron"
    if ch == 1
        unitIdx = [repmat(chanswUnits(ch),ChanUnits(ch),1) [1:ChanUnits(ch)]'];
    elseif ch > 1
        unitIdx = [unitIdx; repmat(chanswUnits(ch),ChanUnits(ch),1) [1:ChanUnits(ch)]'];
    end
end
numUnits = sum(ChanUnits);


display('Aligning AP data on stimulus and response.');
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
    
    %% different setup this time we want the data structure to look like:
    % data(aS).trial(tt).unit(un).times
    
    % loooping over trials
    for tt = 1:nTrials
        % now looping over units
        for un = 1:numUnits
            % getting unit times for the current channel and unit.
            unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==unitIdx(un,1) & ChanUnitTimestamp(:,2)==unitIdx(un,2),3); % in seconds
            
            %% putting the data in a structure.
            data(aS).trial(tt).unit(un).times...
                = unitTimes(unitTimes>(trialStarts(tt)-pre) & unitTimes<trialStarts(tt)+post)...
                - repmat(trialStarts(tt)-pre,length(sum(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post)),1);
            
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
tmp = double(NS3.Data{2});


%% parsing channel labels
if strcmp(patientID,'CUCX2')
    LFPlabels = [repmat('G',1,64)' num2str([1:64]')] ; % anterior ECoG contact to the array.
    LFPlabels = 'G56';
end
cCH = size(LFPlabels,1);

for ch = 1:cCH
    %% TODO: figure out a smart way of allocating which LFP contacts you want
    display('Aligning LFP data on stimulus and response.');
    % for loop to save multiple epochs
    for aS = 1
        % which alignment spot
        switch aS
            case 1
                alignName = 'Cue';
                trialStarts =  trigTimes(trigs>=1 & trigs<28);
                nTrials = sum(trigs>=1 & trigs<28);
                pre = 2;
                post = 5;
            case 2
                alignName = 'Response';
                trialStarts =  trigTimes(trigs>=100 & trigs<=105);
                nTrials = sum(trigs>=100 & trigs<=104);
                pre = 3;
                post = 5;
        end
        
        LFPmat{aS} = zeros(((pre+post)*params.Fs)+1,nTrials);
        for tt = 1:nTrials
            %% tensorizing LFP:
            if floor((trialStarts(tt)+post)*params.Fs)<size(tmp,2)
                LFPmat{aS}(:,tt) = tmp(ch,floor((trialStarts(tt)-pre)*params.Fs):floor((trialStarts(tt)+post)*params.Fs));
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
    
    
    %% [20160622] rejecting the union of the bad trials from both channels.
    rejected = rejectTrials.rejectTheseTrials;
    
    %% now actually rejecting bad trials. [{20161104} I don't actually use this info in this script yet]
    notRejected = ones(1,nTrials);
    notRejected(rejected) = false;
    trialType = trialType.*notRejected;
    trialTypeRJCTD = trialType(logical(notRejected));
    
    
    %% Coherence
    for aS2 = 1
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
        
        
        for t2 = 1:nTrials
            %% calculating avraged cohereograms for each condition
            display(sprintf('calculating spike-field coherence across for trial %d of %d',t2,nTrials))
            % coherence calculation for each trialType
            
            
            
            timeYouWantThisToTake = 'short please';
            if strcmp(timeYouWantThisToTake,'way long')
                %             tic
                %             [C0,phi0,~,~,~,t,f,~,phiStd0,Cerr0] ...
                %                 = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==1)),data(aS2).channel(ch).unit(un).trial(trialType==1), movingWin, params);
                %
                %             [C1a,phi1a,~,~,~,t,f,~,phiStd1a,Cerr1a] ...
                %                 = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==2)),data(aS2).channel(ch).unit(un).trial(trialType==2), movingWin, params);
                %
                %             [C1b,phi1b,~,~,~,t,f,~,phiStd1b,Cerr1b] ...
                %                 = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==3)),data(aS2).channel(ch).unit(un).trial(trialType==3), movingWin, params);
                %
                %             [C2,phi2,~,~,~,t,f,~,phiStd2,Cerr2] ...
                %                 = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==4)),data(aS2).channel(ch).unit(un).trial(trialType==4), movingWin, params);
                %
                %             tspec = linspace(t(1)-pre,t(end)-pre,length(t));
                %             toc
                
            elseif strcmp(timeYouWantThisToTake,'short please')
                % remember data looks like::
                % data(aS).trial(tt).unit(un).times
                
                tic
                %             keyboard
                [C(:,:,t2),phi(:,:,t2),~,~,~,t,f] ...
                    = cohgramcpt(repmat(squeeze(LFPmat{aS2}(:,t2)),1,numUnits),data(aS2).trial(t2).unit, movingWin, params);
                
                tspec = linspace(t(1)-pre,t(end)-pre,length(t));
                toc
                
            end
            %                 end
            
        end % looping over trials
        
        
        savedFileName = sprintf('/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/%s/Data/populationCoherence/%sAlignedCoherence_averagedCohgramsOverUnits_session%d_continuousChannel_G%d.mat',patientID,alignName,sessionNum,ch);
        
        
        %% [20160920] saving coherence data over trials.
        saveData=1;
        if saveData
            save(savedFileName,'C','phi','t','f','tspec')    %,'nullDist','nullPhase'
        end
        
        clear C phi t f
        
        %     set(0,'DefaultFigureRenderer','painters')
        %
        %     imAlpha = 1;
        %     % plotting mean cohereograms over trials.
        %     figure(sessionNum)
        %
        %     for panel = 1:4
        %         % plotting easy coherograms
        %         plotmultipleaxes(panel,4,1,0.04,sessionNum)
        %         hold on
        %         imagesc(tspec,f,squeeze(nanmean(C(:,:,trialType==panel),3))','alphadata',imAlpha)
        %         line([0 0],[params.fpass(1) params.fpass(2)],'color','k','linestyle','--','linewidth',2)
        %         line([nanmedian(rt(trialType(1:nTrials)==panel)./1000) nanmedian(rt(trialType(1:nTrials)==panel)./1000)],[params.fpass(1) params.fpass(2)],'color','k','linestyle','-.','linewidth',2)
        %         hold off
        %         axis xy square tight
        %         xlim([-pre+(movingWin(1)/2) post-(movingWin(1)/2)])
        %         set(gca, 'linewidth',2,'fontsize',14)
        %         colormap(jet)
        %         colorbar('NorthOutside')
        %         set(gca, 'linewidth',2,'fontsize',14)
        %         xlabel('time (seconds)','fontsize',12)
        %         ylabel('LFP Frequency (Hz)','fontsize',12)
        %     end
        %
        %     suptitle(sprintf('SFC for %s LFP channel, AP channel %s, unit %d, aligned on %s.'...
        %         ,LFPlabels,deblank(unitChanLabels(ch,:)),un,alignName));
        %
        %     maximize(sessionNum)
        %     pause(5)
        %     saveas(sessionNum,['/home/elliot/Dropbox/Figs/MSIT-ACC-dlPFC/singleTrialCoherence_UtahArray/coherenceAveragedOverSimultaneouslyRecordedNeurons_session' num2str(sessionNum) '.pdf'])
        
        
        
    end% looping over align spots (Stimulus & response)
end

