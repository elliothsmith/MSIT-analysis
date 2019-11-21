function fieldStats = analyzeMSITfieldsU_appeal(patientID,sessionNum,nevFile,tf)
%ANALYZEMSITFIELDSU analyzes local field potentials recorded on blackrock
%   during the MSIT. This Function is for use with ECoG Data.
%
%   [fieldStats] = analyzeMSITfieldsU(patientID,sessionNum,nevFile,tf)
%   calculates spike-field coherograms for MSIT data in nevFile and its
%   associated ns(3) file. analyzeMSITcoherence saves statistics, and
%   plots results.
%
%   the [tf] input argment is a boolean operator. When tf equals 1,
%       analyzeMSITfields will create spectrograms. When tf equals 0,
%       analyzeMSITfields will create spectra between cue & response.
%
%   Two tips for easy usage:
%       1) place the nev file in the same directory as the nsx files.
%       2) make sure that the nev file with sorted units has its original
%           name at the beginning of the file name.
%

% author: EHS 20160712
% VersionControl: https://github.com/elliothsmith/MSIT-analysis


%% loading data from NEV file
display('loading action potential and local field potential data...')
[dataPath, nvName, nvExt] = fileparts(nevFile);

% defining nsFile.
nsFile = fullfile(dataPath,[nvName(1:19) '.ns3']);

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


%% [20160824] removing practice trial triggers and updating number of trials.
% [trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
% nTrials = sum(trigs==90);


%% parsing behavior &  making a vector of conflict types.
trialType = zeros(1,nTrials);
condition = trigs(trigs>=1 & trigs<=27);
% These are the correct codes. Double Checked on 20160216
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)


%% organizing responses & calculating reaction times
cues = trigTimes(trigs>=1 & trigs<28);
responses = double(trigTimes(trigs>=100 & trigs<=104))';
RTs = responses-cues(1:length(responses));

%% parsing channel labels
macroLabels = deblank({NS3.ElectrodesInfo.Label})';
BFlabels = macroLabels;


%% defining Chronux parameters.
movingWin = [1.2 0.020];
params.Fs = 2e3; % sampling frequency for LFP
params.pad = 1; %
params.fpass = [1 50]; % frequency range of interest
params.tapers = [5 9]; % emphasize smoothing
params.trialave = 0; % average over trials {CHANGES BELOW}
params.err = [2 0.01]; % population error bars


%% //conflict colors//
col0 = [183 30 103]./255;
col1a = [246 139 31]./255;
col1b = [0 166 81]./255;
col2 = [82 79 161]./255;


%% rejecting discharges
[~,rejectTrials] = dischargeRejectionToolU(patientID,sessionNum,nevFile,0)
% defning bad trials
shits = rejectTrials;
tvec = zeros(nTrials,1);
tvec(shits) = 1;
goods = find(~tvec);


%% adjusting trial numbers
trialType = trialType(goods);
RTs = RTs(goods);
nTrials = length(goods);


%% [20160610] denoising LFP
denoiseMethod = 'PCA'
switch denoiseMethod
    case 'PCA'
        display('denoising using PCA...')
        dData = remove1stPC(double(NS3.Data{2}));
        display('...done.')
    case 'notch'
        Wo = 60/(params.Fs/2);  BW = Wo/50;
        [b,a] = iirnotch(Wo,BW);
        %        freqz(b,a);
        for c = 1:size(NS3.Data{2},1)
            display(sprintf('applying notch filter to channel %d',c))
            dData(c,:) = filtfilt(b,a,double(NS3.Data{2}(c,:)));
        end
end


%% building LFP tensor and plotting ERPs and spectrograms.
display('Aligning LFP data on stimulus and response.');
% Channel labels
LFPlabels = [repmat('G',1,64)' num2str([1:64]')] ; % anterior ECoG contact to the array.\
LFPlabels = 'G56';
% cue timing
pre = 2;
post = 4;
% initializing LFPmat
LFPmat = zeros(((pre+post)*params.Fs)+1,length(goods),length(LFPlabels),2); % 2 for cue and response aligned.
% for loop to save multiple epochs

alignName = 'Cue';
trialStarts =  trigTimes(trigs>=1 & trigs<28);
%             nTrials = sum(trigs>=1 & trigs<28);

%% time space for LFP samples.
tsec = linspace(-pre,post,size(LFPmat,1));


% looping over trials.
for tt = 1:length(goods)
    
    
    %% constructing LFP tensor.
    trialLFP = dData(56,floor((trialStarts(goods(tt))-pre)*params.Fs):floor((trialStarts(goods(tt))+post)*params.Fs));
    
    
    %% [20160610] calculating single trial spectrograms here IF THEY ARE NEEDED.
    params.trialave = 0;
    if isequal(mod(tt,50),0)
        display(sprintf('calculated spectrograms for for first %d trials out of %d...',tt,nTrials))
    end
    [Stf(:,:,tt),period] = basewaveERP(trialLFP,params.Fs,params.fpass(1),params.fpass(2),6,0);
    tspec = tsec;
    f = period.^-1;
    
    
    % averaging to collect spectra.
    Sf_encodingWindow(:,tt) = nanmean(abs(Stf(:,tspec>0 & tspec<RTs(tt),tt)),2);
    Sf_baseline(:,tt) = nanmean(abs(Stf(:,tspec>-1.5 & tspec<-0.5,tt)),2);
    
    
end % looping over trials


%% calculating and plotting trial averaged spectrograms.
switch tf
    case 1
        spectralType = 'Spectrograms';
        cLim = [0.8 1.2];
        
        figure(1)
        imagesc(tspec,f,squeeze(nanmean(abs(Stf),3))./repmat(nanmean(Sf_baseline,2),1,size(Stf,2)),cLim);
        colorbar
        axis xy tight square
        xlim([-1 3])
        
end


fName = sprintf('%s_session_%d_Channel_%s_LFP_%sAligned',patientID,sessionNum,deblank(LFPlabels),alignName);
print(1,sprintf('~/data/msit_units/Experiment_II__dlPFC/%s/Figs/spectrograms/%s',patientID,fName),'-dpdf')
close(1)


% exporting some data for stats.
fieldStats.Stf = squeeze(nanmean(abs(Stf),3));
fieldStats.tspec = tspec;
fieldStats.f = f;
fieldStats.Sf_encodingWindow = Sf_encodingWindow;



