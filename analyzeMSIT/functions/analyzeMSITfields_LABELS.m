function LFPlabels = analyzeMSITfields_LABELS(patientID,sessionNum,nevFile,tf)
%ANALYZEMSITFIELDS analyzes local field potentials recorded on blackrock
%   during the MSIT.
%
%   [fieldStats] = analyzeMSITfields(patientID,sessionNum,nevFile,tf)
%   calculates spike-field coherograms for MSIT data in nevFile and its
%   associated ns(3) file. analyzeMSITcoherence saves statistics, and
%   plots results.
%
%   the [tf] input argment flags which type of analyses to perform. 
%       When tf equals 1, analyzeMSITfields will create and save 
%       spectrograms. When tf equals 0, analyzeMSITfields will create and 
%       save spectra between cue & response. When tf equals 2, mean spectra
%       across trials will be created and saved. 
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

% 
% %% organizing important task parameters.
% trigs = NEV.Data.SerialDigitalIO.UnparsedData;
% trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
% TimeRes = NEV.MetaTags.TimeRes;
% nTrials = sum(trigs==90);

% 
% %% [20161017] removing practice trial triggers and updating number of trials.
% if strcmp(patientID,'CUBF09')
%     trigs(95) = 104;
%     [~,trigTimes] = removePracticeTriggers(trigs,trigTimes);
% else
%     [trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
% end
% nTrials = sum(trigs==90);
% 
% 
% %% parsing behavior &  making a vector of conflict types.
% trialType = zeros(1,nTrials);
% condition = trigs(trigs>=1 & trigs<=27);
% % These are the correct codes. Double Checked on 20160216
% trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
% trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
% trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
% trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)
% 
% 
% %% organizing responses & calculating reaction times
% cueTimes = trigTimes(trigs>=1 & trigs<28);
% responseTimes = double(trigTimes(trigs>=100 & trigs<=104));
% RTs = responseTimes(1:length(cueTimes))-cueTimes;
% 
% 
% %% parsing channel labels
% macroLabels = deblank({NS3.ElectrodesInfo.Label})';
% numBFs = length(macroLabels)./8;
% % BFlabels =
% 
% for bf = 1:numBFs
%     tmp = char(deblank(macroLabels((bf)*8,:)));
%     BFlabels{bf} = tmp(1:end-1);
% end
% BFlabels


% %% defining Chronux parameters.
% movingWin = [1.2 0.020];
% params.Fs = 2e3; % sampling frequency for LFP
% params.pad = 1; %
% params.fpass = [1 50]; % frequency range of interest
% params.tapers = [5 9]; % emphasize smoothing
% params.trialave = 0; % average over trials {CHANGES BELOW}
% params.err = [2 0.01]; % population error bars
% 
% 
% %% //conflict colors//
% col0 = [183 30 103]./255;
% col1a = [246 139 31]./255;
% col1b = [0 166 81]./255;
% col2 = [82 79 161]./255;


% %% rejecting discharges
% [~,rejectTrials] = dischargeRejectionTool(patientID,sessionNum,nevFile,0)
% save(fullfile(dataPath,sprintf('%s_rejectedTrialsMSIT_session%d.mat',patientID,sessionNum)),'rejectTrials','nevFile')
% % defning bad trials
% shits = rejectTrials;
% tvec = zeros(nTrials,1);
% tvec(shits) = 1;
% goods = find(~tvec);
% 
% 
% %% adjusting trial numbers
% trialType = trialType(goods);
% RTs = RTs(goods);
% nTrials = length(goods);


% %% [20160610] denoising LFP
% denoiseMethod = 'PCA';
% switch denoiseMethod
%     case 'PCA'
%         display('denoising using PCA...')
%         dData = remove1stPC(double(NS3.Data));
%         display('...done.')
%     case 'notch'
%         display('denoising using a ntoch filter...')
%         Wo = 60/(params.Fs/2);  BW = Wo/50;
%         [b,a] = iirnotch(Wo,BW);
%         %        freqz(b,a);
%         for c = 1:size(NS3.Data,1)
%             display(sprintf('applying notch filter to channel %d',c))
%             dData(c,:) = filtfilt(b,a,double(NS3.Data(c,:)));
%         end
%    case 'notchHighPass'
%         display('denoising using a ntoch filter...')
%         Wo = 60/(params.Fs/2);  BW = Wo/50;
%         [b,a] = iirnotch(Wo,BW);
%         [bh,ah] = ellip(4,0.5, 20, 4/(params.Fs/2),'high');
%         %        freqz(b,a);
%         for c = 1:size(NS3.Data,1)
%             display(sprintf('applying notch filter to channel %d',c))
%             preData(c,:) = filtfilt(b,a,double(NS3.Data(c,:)));
%             dData(c,:) = filtfilt(bh,ah,preData(c,:));
%         end
% end


%% building LFP tensor and plotting ERPs and spectrograms.
display('Aligning LFP data on stimulus and response.');
% Channel labels
LFPlabels = {NS3.ElectrodesInfo.Label};
% % cue timing
% pre = 2;
% post = 4;
% % initializing LFPmat
% LFPmat = zeros(((pre+post)*params.Fs)+1,length(goods),length(LFPlabels),2); % 2 for cue and response aligned.
% for loop to save multiple epochs
% for aS = 1
%     % which alignment spot
%     switch aS
%         case 1
%             alignName = 'Cue';
%             trialStarts =  trigTimes(trigs>=1 & trigs<28);
%             %             nTrials = sum(trigs>=1 & trigs<28);
%         case 2
%             alignName = 'Response';
%             trialStarts =  trigTimes(trigs>=100 & trigs<=105);
%             %             trigTimes(trigs>=100 & trigs<104)
%             % response timing
%             pre = 3;
%             post = 2;
%         case 3
%             alignName = 'Feedback';
%             trialStarts =  trigTimes(trigs>=200 & trigs<=205);
%             % response timing
%             pre = 2;
%             post = 3;
%     end
%     
%     %% time space for LFP samples.
%     tsec = linspace(-pre,post,size(LFPmat,1));
%     
%     % looping over channels.
%     for ch = 1:length(LFPlabels)
%         % looping over trials.
%         for tt = 1:length(goods)
%             %% constructing LFP tensor.
%             try
%                 LFPmat(:,tt,ch,aS) = dData(ch,floor((trialStarts(goods(tt))-pre)*params.Fs):floor((trialStarts(goods(tt))+post)*params.Fs));
%             catch
%                 LFPmat(:,tt,ch,aS) = dData(ch,ceil((trialStarts(goods(tt))-pre)*params.Fs):floor((trialStarts(goods(tt))+post)*params.Fs));
%             end
%             
%             
%             %% [20160610] calculating single trial spectrograms here IF THEY ARE NEEDED.
%             params.trialave = 0;
%             if isequal(mod(tt,50),0)
%                 display(sprintf('Channel %s:: calculated spectrograms for for first %d trials out of %d...',deblank(LFPlabels{ch}),tt,nTrials))
%             end
%             [Stf(:,:,tt),period] = basewaveERP(LFPmat(:,tt,ch,aS),params.Fs,params.fpass(1),params.fpass(2),6,0);
%             tspec = tsec;
%             f = period.^-1;
%             
%             % averaging to collect spectra.
%             Sf_encodingWindow(:,tt) = nanmean(abs(Stf(:,tspec>0 & tspec<RTs(tt),tt)),2);
%             Sf_baseline(:,tt) = nanmean(abs(Stf(:,tspec>-1.5 & tspec<-0.5,tt)),2);
%             
%             
%         end % looping over trials
%         
%         %% calculating and plotting trial averaged spectrograms.
%         switch tf
%             case 1
%                 spectralType = 'Spectrograms';
%                 cLim = [0.8 1.2];
%                 
%                 figure(1)
%                 imagesc(tspec,f,squeeze(nanmean(abs(Stf),3))./repmat(nanmean(Sf_baseline,2),1,size(Stf,2)),cLim);
%                 colorbar
%                 axis xy tight square
%                 xlim([-1 3])
%                 
%         end
%         
%         fName = sprintf('%s_session_%d_Channel_%s_LFP_%sAligned',patientID,sessionNum,deblank(LFPlabels{ch}),alignName);
%         print(1,sprintf('/media/user1/data4TB/data/msit_units/Experiment_I__ACC/%s/Figs/spectrograms/%s',patientID,fName),'-dpdf')
%         close(1)

%         fieldStats(ch).Stf = squeeze(nanmean(abs(Stf),3));
%         fieldStats(ch).tspec = tspec;
%         fieldStats(ch).f = f;
%         fieldStats(ch).Sf_baseline = nanmean(Sf_baseline,2);
%         fieldStats(ch).Sf_encodingWindow = nanmean(Sf_encodingWindow,2);
%         
%         
%     end % looping over lfp channels
% end % looping over align spots (Stimulus & response)
