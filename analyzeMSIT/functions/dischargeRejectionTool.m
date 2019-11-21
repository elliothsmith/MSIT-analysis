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
try
    nsFile = fullfile(dataPath,[nvName(1:19) '.ns3']);
catch
    nsFile = fullfile(dataPath,[nvName(1:8) '*.ns3']);
end

% parsing files and loading data.
if strcmp(nvExt,'.nev')
    NEV = openNEV(nevFile,'read');
    if isequal(exist(nsFile,'file'),0)
        [nsFile,nsPath,~] = uigetfile('*.ns*','Select the correct NSx file');
        NS3 = openNSx(fullfile(nsPath,nsFile));
    else
        NS3 = openNSx(nsFile);
    end
elseif strcmp(nvExt,'.mat')
    load(nevFile);
    if isequal(exist(nsFile,'file'),0)
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
Fs = 2e3;


%% [20160824] removing practice trial triggers and updating number of trials.
[trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
nTrials = sum(trigs==90);


%% timing (seconds)
pre = 2;
post = 3;


%% parsing channel labels
macroLabels = deblank({NS3.ElectrodesInfo.Label})';
numBFs = length(macroLabels)./8;
% BFlabels =

for bf = 1:numBFs
    tmp = char(deblank(macroLabels((bf)*8,:)));
    BFlabels{bf} = tmp(1:end-1);
end
BFlabels


%% NS3.Data => dData (rather than de-noising)
dData = double(NS3.Data);


%% building LFP tensor and plotting ERPs and spectrograms.
display('Aligning LFP data on stimulus and response.');
% initializing LFPmat
LFPlabels = {NS3.ElectrodesInfo.Label};
LFPmat = zeros(((pre+post)*Fs)+1,nTrials,length(LFPlabels),2);
% for loop to save multiple epochs
for aS = 1
    %% going to start by just aligning on CUE.
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>=100 & trigs<=105);
    end
    
    for ch = 1:length(LFPlabels)
        
        if plotFlag
            % figure window.
            ah = figure('Color',[0 0 0]);
            hold on
        end
        
        for tt = 1:nTrials
            
            
            %% [20160622] constructing LFP tensor.
            try
                LFPmat(:,tt,ch,aS) = dData(ch,floor((trialStarts(tt)-pre)*Fs):floor((trialStarts(tt)+post)*Fs));
            catch
                LFPmat(:,tt,ch,aS) = dData(ch,ceil((trialStarts(tt)-pre)*Fs):floor((trialStarts(tt)+post)*Fs));
            end
            
            %% [20160622] time vectors
            tsec = linspace(-pre,post,size(LFPmat,1));
            
            if plotFlag
                %% [20160622] plotting LFP from each trial brushing to reject
                plot(tsec,LFPmat(:,tt,ch,aS)+((tt-1)*1000),'color',rgb('white'))
                axis tight off
            end
            
        end
        
%         keyboard 
        
        [rejectTrialsperChannel(ch).rejectTheseTrials,iprange,fence] = outliers(squeeze(range(LFPmat(:,:,ch,aS))));
        rejectTrialsperChannel(ch).channelLabel = deblank(LFPlabels{ch});
        
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
            
            
            %% saving figures.
            try
                display('saving to Elliot"s dropbox. Thanks!')
                dbPath = sprintf('/home/elliot/Dropbox/MSITunits/Experiment_I__ACC/%s/Figs',patientID);
                fName = sprintf('%s/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
                saveas(gcf,fName, 'pdf')
                close(gcf)
            catch
                maximize(gcf)
                if exist(['./' patientID],'dir')
                    try
                        fName = sprintf('%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
                        saveas(gcf,fName, 'pdf')
                        close(gcf)
                    catch
                        mkdir(sprintf('%s/Figs/',patientID))
                        fName = sprintf('%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
                        saveas(gcf,fName, 'pdf')
                        close(gcf)
                    end
                elseif exist('./Figs','dir')
                    fName = sprintf('%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
                    saveas(gcf,fName, 'pdf')
                    close(gcf)
                else
                    fName = sprintf('%s/Figs/%s_session_%d_Channel_%s_RejectedLFPtrials',dbPath,patientID,sessionNum,rejectTrialsperChannel(ch).channelLabel);
                    saveas(gcf,fName, 'pdf')
                    close(gcf)
                end
            end
            
            display(sprintf('figure saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',fName));
            hold off
            
        end
    end % looping over lfp channels
end % looping over align spots (Stimulus & response)

trialsAcrossChannels = unique(trialsAcrossChannels);
