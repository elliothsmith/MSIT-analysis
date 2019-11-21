function [coherenceStats] = analyzeMSITcoherenceU(patientID,sessionNum,nevFile)
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
params.trialave = 0; % average over trials {CHANGES BELOW}
params.err = [1 0.01]; % population error bars


%% creating neural timing variable
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];
ChanUnitTimestamp(ChanUnitTimestamp(:,2)==255,:) = [];


%% parsing channel labels
chanswUnits = unique(ChanUnitTimestamp(:,1));
nChans = length(chanswUnits);
unitChanLabels = [NEV.ElectrodesInfo.ElectrodeLabel]';

aSvec = 1:2;
display('Aligning AP data on stimulus and response.');
% for loop to save multiple epochs
for aS = aSvec
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
                
                %%~~~~~~~~~~~~~~~~~~~~~~~~FROM O)LD CODE~~~~~~~~~~~~~~~~~~~
                % data(ch).channel(un).unit(tt).times = unitTimes(unitTimes>trialStart-pre & unitTimes<trialStart+post) - repmat(trialStart,length(unitTimes(unitTimes>trialStart-pre & unitTimes<trialStart+post)),1);
                %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                
                %% TODO: Do I need to make sure that all of the times are positive?
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
tmp = double(NS3.Data{2});


%% parsing channel labels
if strcmp(patientID,'CUCX2')
    LFPlabels = 'G56'; % anterior ECoG contact to the array.
end
cCH = size(LFPlabels,1);


%% TODO: figure out a smart way of allocating which LFP contacts you want
display('Aligning LFP data on stimulus and response.');
% for loop to save multiple epochs
for aS = aSvec
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            pre = 2;
            post = 5;
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>=100 & trigs<=105);
            pre = 3;
            post = 5;
    end
    
    LFPmat{aS} = zeros(((pre+post)*params.Fs)+1,nTrials);
    
    
    for tt = 1:nTrials
        
        %%~~~~~~~~~~~~~~~~~~~FROM OLD (WORKING) CODE
        % LFPmat(ch,:,trl) = ECoG2kp(ch,ceil(trialStart(trl)*Fs + win(1)*Fs):ceil(trialStart(trl)*Fs + win(2)*Fs));
        %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        %% tensorizing LFP:
        if floor((trialStarts(tt)+post)*params.Fs)<size(tmp,2)
            LFPmat{aS}(:,tt) = tmp(56,floor((trialStarts(tt)-pre)*params.Fs):floor((trialStarts(tt)+post)*params.Fs));
        end
        
    end
end


%% setting up codes for PSTHs over conflict types.
% These are the correct codes. Double Checked on 20160216
condition = trigs(trigs>=1 & trigs<=27);
trialType = zeros(1,nTrials);
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)


%% Rasters and PSTHs
for aS2 = 2
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
            
            %% excluding timed out trials
            condition = condition(~timeOutTrials);
            trialType = trialType(~timeOutTrials);
            condition = condition(1:nTrials);
            trialType = trialType(1:nTrials);
            
    end
    
    for ch = 1:size(data(aS2).channel,2)
        % [20160607] making sure that spike field coherence is calculated
        % between a spike and teh grid channel G55, which was the channel
        % just anterior to the Utah array.
        tmp = deblank(data(aS2).channel(ch).label);
        
        for un = 1:size(data(aS2).channel(ch).unit,2)
                        
            %                 coherenceDataFile = '';
            %                 [nullDist,nullPhase] = generateCoherenceNullDist(500,squeeze(LFPmat(:,:,cCH,aS2)),data(aS2).channel(ch).unit(un).trial,movingWin,params);
            
            
            %% [20160622] rejecting the union of the bad trials from both channels.
            rejected = rejectTrials.rejectTheseTrials;
            
            
            %% now actually rejecting bad trials.
            notRejected = ones(1,nTrials);
            notRejected(rejected) = false;
%             tmpType = trialType.*notRejected;
            trialTypeRJCTD = trialType(logical(notRejected));
            
            
            if isequal(params.trialave,0)
                %% calculating coherogram for each trial
                %                 display(sprintf('calculating spike-field coherence within trials for\n    LFP channel %s,\n    AP channel %s,\n    unit %d,\n    aligned on %s.'...
                %                     ,strtrim(LFPlabels{cCH}),deblank(unitChanLabels(ch,1:5)),un,alignName));
                
                % caluclate coherence without averaging over trials.
                tic
                [C_alltrials,phi_alltrials,~,~,~,t,f] = cohgramcpt(squeeze(LFPmat{aS2}(:,:)),data(aS2).channel(ch).unit(un).trial, movingWin, params);
                tspec = linspace(-pre,post,length(t));
                toc
                
                savedFileName = sprintf('/media/user1/data4TB/msit_units/Experiment_II__dlPFC/%s/Data/coherence/%sAlignedCoherence_singleTrialCohgrams_session%d_unit%d_channel%d_continuousChannel%s.mat',patientID,alignName,sessionNum,un,ch,LFPlabels);
                save(savedFileName,'C_alltrials','phi_alltrials','t','f','tspec',...
                    'pre','post','-v7.3')
                
%                 [nullDist,nullPhase] = generateCoherenceNullDist(500, squeeze(LFPmat{aS2}(:,:)), data(aS2).channel(ch).unit(un).trial, movingWin, params);
                
                % plotting mean cohereograms over trials.
                figure(ch*1000+un)
                
                % plotting easy coherograms
                plotmultipleaxes(1,4,1,0.04,ch*1000+un)
                imagesc(tspec,f,squeeze(nanmean(C_alltrials(:,:,trialType==1),3))',[minval maxval]);
                line([0 0],[params.fpass(1) params.fpass(2)],'color','k','linestyle','--','linewidth',2)
                line([median(rt(trialType(1:nTrials)==1)./1000) median(rt(trialType(1:nTrials)==1)./1000)],[params.fpass(1) params.fpass(2)],'color','k','linestyle','-.','linewidth',2)
                axis xy square tight
                xlim([-pre+(movingWin(1)/2) post-(movingWin(1)/2)])
                set(gca, 'linewidth',2,'fontsize',14)
                colormap(jet)
                set(gca, 'linewidth',2,'fontsize',14)
                xlabel('time (seconds)','fontsize',12)
                ylabel('LFP Frequency (Hz)','fontsize',12)
                title('none')
                
                % plotting easy coherograms
                plotmultipleaxes(2,4,1,0.04,ch*1000+un)
                imagesc(tspec,f,squeeze(nanmean(C_alltrials(:,:,trialType==2),3))',[minval maxval]);
                line([0 0],[params.fpass(1) params.fpass(2)],'color','k','linestyle','--','linewidth',2)
                line([median(rt(trialType(1:nTrials)==2)./1000) median(rt(trialType(1:nTrials)==2)./1000)],[params.fpass(1) params.fpass(2)],'color','k','linestyle','-.','linewidth',2)
                axis xy square tight
                xlim([-pre+(movingWin(1)/2) post-(movingWin(1)/2)])
                set(gca, 'linewidth',2,'fontsize',14)
                colormap(jet)
                set(gca, 'linewidth',2,'fontsize',14)
                xlabel('time (seconds)','fontsize',12)
                ylabel('LFP Frequency (Hz)','fontsize',12)
                title('spatial')
                
                % plotting easy coherograms
                plotmultipleaxes(3,4,1,0.04,ch*1000+un)
                imagesc(tspec,f,squeeze(nanmean(C_alltrials(:,:,trialType==3),3))',[minval maxval]);
                line([0 0],[params.fpass(1) params.fpass(2)],'color','k','linestyle','--','linewidth',2)
                line([median(rt(trialType(1:nTrials)==3)./1000) median(rt(trialType(1:nTrials)==3)./1000)],[params.fpass(1) params.fpass(2)],'color','k','linestyle','-.','linewidth',2)
                axis xy square tight
                xlim([-pre+(movingWin(1)/2) post-(movingWin(1)/2)])
                set(gca, 'linewidth',2,'fontsize',14)
                colormap(jet)
                set(gca, 'linewidth',2,'fontsize',14)
                xlabel('time (seconds)','fontsize',12)
                ylabel('LFP Frequency (Hz)','fontsize',12)
                title('distractor')
                
                % plotting hard coherograms
                plotmultipleaxes(4,4,1,0.04,ch*1000+un)
                imagesc(tspec,f,squeeze(nanmean(C_alltrials(:,:,trialType==4),3))',[minval maxval]);
                line([0 0],[params.fpass(1) params.fpass(2)],'color','k','linestyle','--','linewidth',2)
                line([median(rt(trialType(1:nTrials)==4)./1000) median(rt(trialType(1:nTrials)==4)./1000)],[params.fpass(1) params.fpass(2)],'color','k','linestyle','-.','linewidth',2)
                axis xy square tight
                xlim([-pre+(movingWin(1)/2) post-(movingWin(1)/2)])
                set(gca, 'linewidth',2,'fontsize',14)
                colormap(jet)
                colorbar('NorthOutside')
                set(gca, 'linewidth',2,'fontsize',14)
                xlabel('time (seconds)','fontsize',12)
                ylabel('LFP Frequency (Hz)','fontsize',12)
                title('both')
                                
                
            elseif isequal(params.trialave,1)
                %% calculating avraged cohereograms for each condition
                display(sprintf('calculating spike-field coherence across trials for\n    LFP channel %s,\n    AP channel %s,\n    unit %d,\n    aligned on %s.'...
                    ,LFPlabels,strtrim(unitChanLabels(ch,1:5)),un,alignName));
                % coherence calculation for each trialType
                
                savedFileName = sprintf('/home/elliot/data/msit_units/%s/Data/coherence/%sAlignedCoherence_averagedCohgrams_session%d_unit%d_channel%d_continuousChannel%s.mat',patientID,alignName,sessionNum,un,ch,LFPlabels);
                
                
                timeYouWantThisToTake = 'way long';
                if strcmp(timeYouWantThisToTake,'way long')
                    tic
                    [C0,phi0,~,~,~,t,f,~,phiStd0,Cerr0] ...
                        = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==1)),data(aS2).channel(ch).unit(un).trial(trialType==1), movingWin, params);
                    
                    [C1a,phi1a,~,~,~,t,f,~,phiStd1a,Cerr1a] ...
                        = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==2)),data(aS2).channel(ch).unit(un).trial(trialType==2), movingWin, params);
                    
                    [C1b,phi1b,~,~,~,t,f,~,phiStd1b,Cerr1b] ...
                        = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==3)),data(aS2).channel(ch).unit(un).trial(trialType==3), movingWin, params);
                    
                    [C2,phi2,~,~,~,t,f,~,phiStd2,Cerr2] ...
                        = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==4)),data(aS2).channel(ch).unit(un).trial(trialType==4), movingWin, params);
                    
                    tspec = linspace(t(1)-pre,t(end)-pre,length(t));
                    toc
                    
                elseif strcmp(timeYouWantThisToTake,'short please')
                    tic
                    [C0,phi0,~,~,~,t,f] ...
                        = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==1)),data(aS2).channel(ch).unit(un).trial(trialType==1), movingWin, params);
                    
                    [C1a,phi1a,~,~,~,t,f] ...
                        = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==2)),data(aS2).channel(ch).unit(un).trial(trialType==2), movingWin, params);
                    
                    [C1b,phi1b,~,~,~,t,f] ...
                        = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==3)),data(aS2).channel(ch).unit(un).trial(trialType==3), movingWin, params);
                    
                    [C2,phi2,~,~,~,t,f] ...
                        = cohgramcpt(squeeze(LFPmat{aS2}(:,trialType==4)),data(aS2).channel(ch).unit(un).trial(trialType==4), movingWin, params);
                    
                    tspec = linspace(t(1)-pre,t(end)-pre,length(t));
                    toc
                    
                end
                %                 end
                
                %% [20160920] saving coherence data over trials.
                saveData=1;
                if saveData
                    save(savedFileName,'C0','C1a','C1b','C2',...
                        'phi0','phi1a','phi1b','phi2',...
                        'phiStd0','phiStd1a','phiStd1b','phiStd2',...
                        'Cerr0','Cerr1a','Cerr1b','Cerr2','t','f',...
                        'tspec','-v7.3')    %,'nullDist','nullPhase'
                    
                    statAlpha = median([min(Cerr0) min(Cerr1a) min(Cerr1b) min(Cerr2)]);
                else
                    statAlpha = 0.2;
                    
                    Cerr0 = 0.05;
                    Cerr1a = 0.05 ;
                    Cerr1b = 0.05;
                    Cerr2 = 0.05;
                end
                
                maxval = max([max(max(C0)) max(max(C1a)) max(max(C1b)) max(max(C2))]);
                imAlpha = 0.3;
            
        end % looping over units for each channel and align spot
    end % looping over channels for each align spot
end % looping over align spots (Stimulus & response)


end % end of function

