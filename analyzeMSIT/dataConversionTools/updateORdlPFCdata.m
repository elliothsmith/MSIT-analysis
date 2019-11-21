%% [0628] converts dlPFC data to a NEV format that will allow for similar analyses to the 
dataDir = '~/Dropbox/MSITunits/OR';
preprocessedDir = '/Users/elliotsmith/Dropbox/MSITunits';

dirlist = dir(dataDir);

% looping over patients
for fl = 4:length(dirlist)
    patientID = dirlist(fl).name(1:3);
    
    display(sprintf('Converting data for patient %s...',patientID))
            
    % entering patient-specific directories to get AP data.
    cd([dataDir '/' dirlist(fl).name])
    
    % grabbing trial markers.
    load([patientID '_MSIT_OR_trialmarkers.mat'])

    
    %% loading neural data
    APchannels = dir([patientID '-Channel*.csv']);
    ChanUnitTimestamp = [];
    for ch = 1:length(APchannels)
            spikemat = importdata(APchannels(ch).name);
            if isstruct(spikemat)
                tmp = spikemat.data(:,1:3);
                ChanUnitTimestamp = cat(1,ChanUnitTimestamp,tmp);
            else
                ChanUnitTimestamp = spikemat(:,1:3);
            end
    end
    
    % going back to the BHVs to get response and feedback data.
    switch patientID
        case {'C10'}
            % load file for c17
            load('behav_one.mat')
            startTrial = 1;
        case {'C17'}
            % load file for c17
            load('20140821_run1.mat')
            behav_one = BHV;
            startTrial = 1;
        case {'C18'}
            % load file for c18
            load('20140911_run1.mat')
            behav_one = BHV;
            startTrial = 23;
        case {'C20'}
            % load file for c20
            load('c20_ORMSIT_BHV.mat')
            behav_one = BHV;
            startTrial = 10; % This is the trial the recording started on
        case {'C21'}
            % load file for c21
            load('c21_ORMSIT_BHV.mat')
            behav_one = BHV;
            startTrial = 1;
        case {'C23'}
            % load file for c21
            load('c23_bhv.mat')
            behav_one = BHV;
            startTrial = 21;
        case {'C26'}
            % load file for c27
            load('c26_bhv.mat')
            behav_one = bhv;
            startTrial = 57;
        case {'C30'}
            % load file for c30
            load('C30_MSIT_BHV.mat')
            behav_one = bhv;
            startTrial = 1;
        case {'C33'}
            % load file for c30
            load('C33_MSIT_OR_LFP_BHV.mat')
            behav_one = bhv;
            startTrial = 1;
    end
    
    
    %% going ahead and pulling out the full trial structure from the bhv. 
    %% note to self: Never ever use monkeylogic again...
    nTrials = length(trialMarkers.responses);
    if strcmp(patientID,'C17')
        MLtrialIdcs = startTrial:startTrial+nTrials-2;
        nTrials = length(trialMarkers.responses)-1;
    else
        MLtrialIdcs = startTrial:startTrial+nTrials-1;
    end
        
    % putting monkeylogic codes in a matrix. 
    codes = zeros(16,nTrials);
    for tt = 1:nTrials
        try
            if length(behav_one.CodeNumbers{tt})<16
                codes(1:length(behav_one.CodeNumbers{tt}),tt) = behav_one.CodeNumbers{tt};
            else
                codes(:,tt) = behav_one.CodeNumbers{tt};
            end
        catch
            display(sprintf('patient %s had %d more trials than there is behavioral data for',patientID,nTrials-length(behav_one.CodeNumbers)));
        end
    end
    
    % arranging marker times.
    A = cat(1,trialMarkers.trialStarts-500, trialMarkers.trialStarts, ...
        trialMarkers.fixations, trialMarkers.cues, trialMarkers.responses, ...
        trialMarkers.responses);
    markerTimes = reshape(A,size(A,1)*size(A,2),1);
    
    % arranging marker values.
    startVals0 = repmat(90,1,nTrials);
    startVals1 = repmat(91,1,nTrials);
    fixVals = repmat(92,1,nTrials);
    cueVals = behav_one.ConditionNumber(MLtrialIdcs)';
    
    
    %% response processing. 
    respVals = codes(9,:)+repmat(43,1,nTrials);
    
    
    %% error and feedback processing
    incorrs = behav_one.TrialError(MLtrialIdcs)==6;
    toosoon = behav_one.TrialError(MLtrialIdcs)==5;
    timedout = behav_one.TrialError(MLtrialIdcs)==1 | behav_one.TrialError(MLtrialIdcs)==2;
    fbPresent = logical(codes(5,:)-repmat(48,1,nTrials));
    
    fbVals = repmat(204,1,nTrials);
    fbVals(fbPresent) = 200;    
    for fb = 1:length(fbVals)
        % for error trials, +1::
        if ~isequal(incorrs(tt),0) & fbVals==204
            fbVals(fb) = 205;
        elseif ~isequal(incorrs(tt),0) & fbVals==200
            fbVals(fb) = 205;
        end
    end
    fbVals(toosoon) = 203;    
    fbVals(timedout) = 202;    

    % putting all the values in a vector. 
    B = cat(1, startVals0, startVals1, fixVals, cueVals(1:nTrials), respVals, fbVals);
    markerVals = reshape(B,size(B,1)*size(B,2),1);
    
    
    %% saving NEVs
    NEV.MSITextras.PatientID = patientID;
    NEV.MSITextras.codes = codes;
    NEV.MSITextras.behavioralExtras = behav_one;
    NEV.Data.SerialDigitalIO.UnparsedData = markerVals;  
    NEV.Data.Spikes.Electrode = ChanUnitTimestamp(:,1);
    NEV.Data.Spikes.Unit = ChanUnitTimestamp(:,2);
    NEV.Data.Spikes.TimeStamp = ChanUnitTimestamp(:,3);
    NEV.Data.SerialDigitalIO.UnparsedData = markerVals;
    if strcmp(patientID,'C17')
        NEV.Data.SerialDigitalIO.TimeStampMS = markerTimes(1:length(markerVals));
    else
        NEV.Data.SerialDigitalIO.TimeStampMS = markerTimes;
    end
    NEV.Data.SerialDigitalIO.TimeStampSec = markerTimes./1000;

    save([dataDir '/' patientID '/' dirlist(fl).name '_reconstructedNEV.mat'],'NEV','-v7.3')
    
    clearvars -except dataDir preprocessedDir dirlist fl
    
    
end
