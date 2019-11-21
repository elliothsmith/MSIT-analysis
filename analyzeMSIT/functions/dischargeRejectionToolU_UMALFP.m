function [rejectTrialsperChannel,trialsAcrossChannels] = dischargeRejectionTool(patientID,sessionNum,nevFile,plotFlag)
%DISCHARGREJECTIONTOOL Tool for rejecting trials with epileptiform discharges.
%
%   [rejectTrialsperChannel,trialsAcrossChannels] = dischargeRejectionTool
%       (patientID,sessionNum,nevFile) finds outliers based on their
%       maximum and mimimum values.
%
%   optional input argument: plotFlag == 0 will not plot rejected trials.
%
%   Two tips for easy usage:
%       1) place the nev file in the same directory as the nsx files.
%       2) make sure that the nev file with sorted units has its original
%           name at the beginning of the file name.
%


% author: ElliotHSmith (https://github.com/elliothsmith/MSIT-analysis)

% default argument for plotting.
if nargin==3
    plotFlag = 1;
end


%% loading data from NEV file
display('loading action potential and local field potential data...')
[dataPath, nvName, nvExt] = fileparts(nevFile);

% defining nsFile.
nsFile = fullfile(dataPath,[nvName(1:19) '.ns3']);
if ~exist(nsFile,'file')
    nsList = dir(fullfile(dataPath,[nvName(1:8) '*.ns3']));
    nsFile = fullfile(dataPath,nsList(1).name);
end

% parsing files and loading AP data.
if strcmp(nvExt,'.nev')
    NEV = openNEV(nevFile,'read');
elseif strcmp(nvExt,'.mat')
    load(nevFile);
end

% loading LFP data.
if ~exist(nsFile,'file')
    [nsFile,nsPath,~] = uigetfile('*.ns*','Select the correct NSx file');
    %% loading LFP
    load([nsFile(1:end-4) '_downsampledMeanUMALFP.mat'])
else
    load([nsFile(1:end-4) '_downsampledMeanUMALFP.mat'])
end


%% updatnig variable name for UMA data.
dData = dsData;
clear dsData


%% organizing important task parameters.
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);
Fs = 500;


%% [20160824] removing practice trial triggers and updating number of trials.
[trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
nTrials = sum(trigs==90);


%% timing (seconds)
pre = 2;
post = 3;


%% parsing channel labels
if strcmp(patientID,'CUCX2')
    LFPlabels = 'UMALFP'; % anterior ECoG contact to the array.
end


%% building LFP tensor and plotting ERPs and spectrograms.
display('Aligning LFP data on stimulus and response.');
% initializing LFPmat
LFPmat = zeros(((pre+post)*Fs)+1,nTrials,length(LFPlabels),2);
% for loop to save multiple epochs
for aS = 1
    %% going to start by just aligning on CUE.
    alignName = 'Cue';
    trialStarts =  trigTimes(trigs>=1 & trigs<28);
    nTrials = sum(trigs>=1 & trigs<28);
    
    
    for ch = 1 % the contact just anterior to the array
        
        if plotFlag
            % figure window.
            ah = figure('Color',[0 0 0]);
            hold on
        end
        
        for tt = 1:nTrials
            
            if floor((trialStarts(tt)+post)*Fs)<length(dData)
                
                if floor((trialStarts(tt)+post)*Fs)<size(LFPmat,1) || floor((trialStarts(tt)-pre)*Fs)>0
                    LFPmat(:,tt,aS) = dData(floor((trialStarts(tt)-pre)*Fs):floor((trialStarts(tt)+post)*Fs));
                end
                
                
                %% [20160622] time vectors
                tsec = linspace(-pre,post,size(LFPmat,1));
                
                if plotFlag
                    %% [20160622] plotting LFP from each trial brushing to reject
                    plot(tsec,LFPmat(:,tt,ch,aS)+((tt-1)*1000),'color',rgb('white'))
                    axis tight off
                end
            else
                nTrials = tt-1;
            end
        end
        rejectTrialsperChannel(ch).rejectTheseTrials = outliers(squeeze(range(LFPmat(:,:,ch,aS))));
        rejectTrialsperChannel(ch).channelLabel = LFPlabels; % may need to edit this for many channels if we do another grid case.
        
        % saving all of the bad trials
        if isequal(ch,1)
            trialsAcrossChannels = outliers(squeeze(range(LFPmat(:,:,ch,aS))));
        else
            trialsAcrossChannels = cat(2,trialsAcrossChannels,outliers(squeeze(range(LFPmat(:,:,ch,aS)))));
        end
        
        %%
        if plotFlag
            for rjct = rejectTrialsperChannel(ch).rejectTheseTrials
                %% [20160622] plotting LFP from each trial brushing to reject
                plot(tsec,LFPmat(:,rjct,ch,aS)+((rjct-1)*1000),'color',rgb('springgreen'))
            end
            maximize(gcf)
            title(sprintf('Channel: %s. || Rejected Trials in Green',rejectTrialsperChannel(ch).channelLabel),'color',rgb('springgreen'))
            
            
            %             %% saving figures.
            %             try
            %                 display('saving to Elliot"s dropbox. Thanks!')
            %                 dbPath = sprintf('/home/elliot/Dropbox/MSITunits_emu/%s/Figs',patientID);
            %                 fName = sprintf('%s/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
            %                 saveas(gcf,fName, 'pdf')
            %                 close(gcf)
            %             catch
            %                 maximize(gcf)
            %                 if exist(['../../' patientID],'dir')
            %                     try
            %                         fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
            %                         saveas(gcf,fName, 'pdf')
            %                         close(gcf)
            %                     catch
            %                         mkdir(sprintf('../../%s/Figs/',patientID))
            %                         fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
            %                         saveas(gcf,fName, 'pdf')
            %                         close(gcf)
            %                     end
            %                 elseif exist('../../Figs','dir')
            %                     fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
            %                     saveas(gcf,fName, 'pdf')
            %                     close(gcf)
            %                 else
            %                     fName = sprintf('../../%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
            %                     saveas(gcf,fName, 'pdf')
            %                     close(gcf)
            %                 end
            %             end
            %
            %             display(sprintf('figure saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',fName));
            hold off
        end
    end % looping over lfp channels
end % looping over align spots (Stimulus & response)

trialsAcrossChannels = unique(trialsAcrossChannels);
