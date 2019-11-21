function [aptep] = analyzeMSITaptepU(patientID,sessionNum,nevFile)
%ANALYZEMSITAPTEPU creates action potential-triggered evoked potentials.
%
%   [aptepStats] = analyzeMSITaptepU(patientID,sessionNum,nevFile)
%   calculates spike-field coherograms for MSIT data in nevFile and its
%   associated ns(3) file. analyzeMSITcoherence saves statistics, and
%   plots results.
%
%   Two tips for easy usage:
%       1) place the nev file in the same directory as the nsx files.
%       2) make sure that the nev file with sorted units has its original
%           name at the beginning of the file name.
%           (e.g. ./yyyymmdd-hhmmss-001-anythingYouWantHere.nev)
%


% author: ElliotHSmith (https://github.com/elliothsmith/MSIT-analysis)


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


%% parsing behavior &  making a vector of conflict types.
trialType = zeros(1,nTrials);
condition = trigs(trigs>=1 & trigs<=27);
% These are the correct codes. Double Checked on 20160216
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)


%% organizing responses & calculating reaction times
cueTimes = trigTimes(trigs>=1 & trigs<28);
responseTimes = double(trigTimes(trigs>=100 & trigs<=104));
RTs = responseTimes(1:length(cueTimes))-cueTimes;


%% timing (seconds)
pre = 3;
post = 5;
spikeWin = 3; % seconds on either side of action potential.


%% creating neural timing variable
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];
ChanUnitTimestamp(ChanUnitTimestamp(:,2)==255,:) = [];


%% parsing micro channel labels
chanswUnits = unique(ChanUnitTimestamp(:,1));
nChans = length(chanswUnits);
unitChanLabels = [NEV.ElectrodesInfo.ElectrodeLabel]';


%% defining Chronux parameters.
movingWin = [2 0.020];
params.Fs = 2e3; % sampling frequency for LFP
params.pad = 1; %
params.fpass = [0 25]; % frequency range of interest
params.tapers = [5 9]; % emphasize smoothing for the spikes
params.trialave = 1; % average over trials
params.err = [2 0.01]; % population error bars


%% //conflict colors//
col0 = [183 30 103]./255;
col1a = [246 139 31]./255;
col1b = [0 166 81]./255;
col2 = [82 79 161]./255;


%% rejecting discharges
[~,rejectTrials] = dischargeRejectionToolU(patientID,sessionNum,nevFile,0)
save(fullfile(dataPath,sprintf('%s_rejectedTrialsMSIT_session%d.mat',patientID,sessionNum)),'rejectTrials','nevFile')
% defning bad trials
shits = rejectTrials;
tvec = zeros(nTrials,1);
tvec(shits) = 1;
goods = find(~tvec);


%% adjusting trial numbers
trialType = trialType(goods);
RTs = RTs(goods);
nTrials = length(goods);
if sessionNum>2; nTrials = nTrials-1; trialType = trialType(1:end-1); end


%%  alinging LFP data on the same time base as AP data.
tmp = double(NS3.Data{2});


%% parsing channel labels
if strcmp(patientID,'CUCX2')
    LFPlabels = 'G56'; % anterior ECoG contact to the array.
end
cCH = 56;


%% [20160610] denoising LFP
denoiseMethod = 'notch';
switch denoiseMethod
%     case 'PCA'
%         display('denoising using PCA...')
%         dData = remove1stPC(double(NS3.Data));
%         display('...done.')
    case 'notch'
        display('denoising using a ntoch filter...')
        Wo = 60/(params.Fs/2);  BW = Wo/50;
        [b,a] = iirnotch(Wo,BW);
        %        freqz(b,a);
        display(sprintf('applying notch filter to channel %s',LFPlabels))
        dData = filtfilt(b,a,double(NS3.Data{2}(cCH,:)));
%     case 'notchHighPass'
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
end


%% main loop.
display('Aligning AP data on stimulus and response.');
% for loop to save multiple epochs
for aS = 1
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>=100 & trigs<=105);
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
                % updating user.
                if mod(tt,50)==0
                    display(sprintf('built LFP tensor aligned on AP times for trial %d of %d for unit %d of %d for channel %d of %d.',tt,nTrials,un,nUnits,ch,nChans));
                end
                
                % [20160615] aligning data on each cue without subtracting
                % the time of the cue. => just making a spike time struct.
                data(aS).channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>(trialStarts(goods(tt))-pre) & unitTimes<trialStarts(goods(tt))+post);
                
                % this time I don't want to align on the cue onset, but I want
                % the spikes that are just after the cue.
                % - repmat(trialStarts(tt)-pre,length(sum(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post)),1);
                whichSpikes = [0 3];
                
                
                %% defining the AP times on which to align.
                relevantEvents{tt} = data(aS).channel(ch).unit(un).trial(tt).times(data(aS).channel(ch).unit(un).trial(tt).times > trialStarts(goods(tt))+whichSpikes(1) & data(aS).channel(ch).unit(un).trial(tt).times <= trialStarts(goods(tt))+whichSpikes(2));
                
            end % for trials.
            
            
            % building spike matrix.
            for tt2 = 1:nTrials
                if isempty(relevantEvents{tt2})
                    tmpAPmat(:,tt2) = nan(2*spikeWin*params.Fs+1,1);
                else
                    for zz = 1:length(relevantEvents{tt2})
                        try
                            tmptmpAPmat(:,zz) = dData(1,round(relevantEvents{tt2}(zz)-spikeWin)*params.Fs:round(relevantEvents{tt2}(zz)+spikeWin)*params.Fs);
                        end
                    end
                    tmpAPmat(:,tt2) = mean(tmptmpAPmat(:,zz),2);
                    clear tmptmpAPmat
                end
            end
            
            
            %% mean apteps for specific categories of trials.
            aptepAll = tmpAPmat;
            tsec = linspace(-spikeWin,spikeWin,size(aptepAll,1));
            
            
            %% [20161121] putting aptep data in a struct for saving.
            %                 keyboard
            aptep.data{ch,un} = tmpAPmat;
            aptep.tsec{ch,un} = tsec;
            aptep.alignName{ch,un} = alignName;
            aptep.spikeChannelLabel{ch,un} = data(aS).channel(ch).label;
            aptep.spikeChannel{ch,un} = ch;
            aptep.spikeunit{ch,un} = un;
            aptep.lfpChannelLlabel{ch,un} = 'G56';
            aptep.trialLabels{ch,un} = trialType;
            
            
            %% plotting the mean aptep
            set(0, 'DefaultFigureRenderer', 'painters');
            figure (ch*1000+un)
            
            % plotting apteps
            plotmultipleaxes(1,1,3,0.05,ch*1000+un)
            hold on
            plot(tsec,squeeze(nanmean(tmpAPmat,2)),'color','k','linewidth',1)
            plot(tsec,squeeze(nanmean(tmpAPmat(:,trialType==2),2)),'color',col1a,'linewidth',1)
            plot(tsec,squeeze(nanmean(tmpAPmat(:,trialType==3),2)),'color',col1b,'linewidth',1)
            plot(tsec,squeeze(nanmean(tmpAPmat(:,trialType==1),2)),'color',col0,'linewidth',1)
            plot(tsec,squeeze(nanmean(tmpAPmat(:,trialType==4),2)),'color',col2,'linewidth',1)
            text(2,20,'-:medial, --:lateral','fontsize',16,'fontweight','bold')
            hold off
            axis tight
            xlim([-1 1])
            set(gca,'linewidth',1.5,'fontsize',14)
            
            
            %% [20160622] doing stats to find significant regions
            % sliding GLM with 200 ms window and 50 ms step size.
            lmWin = 0.500*params.Fs; % in samples
            lmStep = 0.050*params.Fs;
            lmRange = size(aptepAll,1);
            numSteps = floor(lmRange/lmStep)-round(lmWin/lmStep);
            
            for stp = 1:numSteps
                stpIdx = ((stp-1)*lmStep)+1:((stp-1)*lmStep)+lmWin;
                [p,aptepAnova(stp).anovatab,aptepAnova(stp).stats] = anova1(aptepAll(stpIdx,:),trialType','off');
                
                % plotting significance.
                if p<0.01
                    patch([tsec(stpIdx(1)) tsec(stpIdx(end)) tsec(stpIdx(end)) tsec(stpIdx(1))],[19 19 20 20],rgb('dimgray'),'edgecolor','none')
                end
                
            end
            
            
            %% calculating spectral decompositions of APtEPs
%             display('now calculating spectra for the same unit and LFP channel...')
%             [S0,t,f] = mtspecgramc(tmpAPmat(:,trialType==1),movingWin,params);
%             [S1a,t,f] = mtspecgramc(tmpAPmat(:,trialType==2),movingWin,params);
%             [S1b,t,f] = mtspecgramc(tmpAPmat(:,trialType==3),movingWin,params);
%             [S2,t,f] = mtspecgramc(tmpAPmat(:,trialType==4),movingWin,params);
%             tspec = linspace(-spikeWin,spikeWin,length(t));
%             %                 keyboard
            
            
%             %% plotting spectral decomposition of APtEPs.
%             plotmultipleaxes(2,4,2,0.05,ch*1000+un)
%             hold on
%             imagesc(tspec,f,real(S0./repmat(f.^-2,size(S0,1),1))')
%             line([0 0],params.fpass,'color','k','linestyle','--','linewidth',2)
%             text(-1,params.fpass(2)-10,sprintf('none%s','G56'))
%             hold off
%             axis tight xy square
%             set(gca,'linewidth',1.5,'fontsize',14)
%             plotmultipleaxes(4,4,2,0.05,ch*1000+un)
%             hold on
%             imagesc(tspec,f,real(S1a./repmat(f.^-2,size(S1a,1),1))')
%             line([0 0],params.fpass,'color','k','linestyle','--','linewidth',2)
%             text(-1,params.fpass(2)-10,sprintf('simon%s','G56'))
%             hold off
%             axis tight xy square
%             set(gca,'linewidth',1.5,'fontsize',14)
%             plotmultipleaxes(6,4,2,0.05,ch*1000+un)
%             hold on
%             imagesc(tspec,f,real(S1b./repmat(f.^-2,size(S1b,1),1))')
%             line([0 0],params.fpass,'color','k','linestyle','--','linewidth',2)
%             text(-1,params.fpass(2)-10,sprintf('eriksen%s','G56'))
%             hold off
%             axis tight xy square
%             set(gca,'linewidth',1.5,'fontsize',14)
%             plotmultipleaxes(8,4,2,0.05,ch*1000+un)
%             hold on
%             imagesc(tspec,f,real(S2./repmat(f.^-2,size(S2,1),1))')
%             line([0 0],params.fpass,'color','k','linestyle','--','linewidth',2)
%             text(-1,params.fpass(2)-10,sprintf('both%s','G56'))
%             hold off
%             axis tight xy square
%             set(gca,'linewidth',1.5,'fontsize',14)
%             
%             colormap jet
%             
            
            %% [20161121] saving aptep data.
            if exist('./Data','dir')
                try
                    sfName = sprintf('./Data/APtEP/%s_session_%d_APtEPdata03_denoisedWith%s.mat',patientID,sessionNum,denoiseMethod);
                catch
                    mkdir('./Data/APtEP/')
                    sfName = sprintf('./Data/APtEP/%s_session_%d_APtEPdata03_denoisedWith%s.mat',patientID,sessionNum,denoiseMethod);
                end
            elseif exist('./Data/APtEP/','dir')
                sfName = sprintf('../../Data/APtEP/%s_session_%d_APtEPdata03_denoisedWith%s.mat',patientID,sessionNum,denoiseMethod);
            else
                sfName = sprintf('%s_session_%d_APtEPdata03_denoisedWith%s.mat',patientID,sessionNum,denoiseMethod);
            end
            
            save(sfName,'aptep','-v7.3')
            display(sprintf('Data saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',sfName));
            
            
            handleF = ch*1000+un;
            
            %         %% saving figures.
            %         try
            %             display('saving to Elliot"s dropbox. Thanks!')
            %             dbPath = sprintf('/home/elliot/Dropbox/MSITunits_emu/%s/Figs/LFP/',patientID);
            %             fName = sprintf('%s/%s_session_%d_Channel_%s_LFP_%sAligned',dbPath,patientID,sessionNum,deblank(LFPlabels{ch}),alignName);
            %             maximize(handleF)
            %             pause(3)
            %             saveas(handleF,fName,'pdf')
            %             close(handleF)
            %         catch
%             maximize(handleF)
%             pause(3)
%             if exist('./Figs','dir')
%                 try
%                     fName = sprintf('./Figs/APtEP/%s_Data03_session_%d_channel%d_unit%d_LFPChannel_%s_APtEP_denoisedWith%s',patientID,sessionNum,ch,un,'G56',denoiseMethod);
%                     %                     saveas(handleF,fName,'pdf')
%                     print(fName,'-dpdf','-fillpage')
%                     close(handleF)
%                 catch
%                     mkdir('./Figs/APtEP/')
%                     fName = sprintf('./Figs/APtEP/%s_Data03_session_%d_channel%d_unit%d_LFPChannel_%s_APtEP_denoisedWith%s',patientID,sessionNum,ch,un,'G56',denoiseMethod);
%                     %                     saveas(handleF,fName,'pdf')
%                     print(fName,'-dpdf','-fillpage')
%                     close(handleF)
%                 end
%             elseif exist('./Figs/APtEP/','dir')
%                 fName = sprintf('../../Figs/APtEP/%s_Data03_session_%d_channel%d_unit%d_LFPChannel_%s_APtEP_denoisedWith%s',patientID,sessionNum,ch,un,'G56',denoiseMethod);
%                 %                     saveas(handleF,fName,'pdf')
%                 print(fName,'-dpdf','-fillpage')
%                 close(handleF)
%             else
%                 fName = sprintf('%s_Data03_session_%d_channel%d_unit%d_LFPChannel_%s_APtEP_denoisedWith%s',patientID,sessionNum,ch,un,'G56',denoiseMethod);
%                 %                     saveas(handleF,fName,'pdf')
%                 print(fName,'-dpdf','-fillpage')
%                 close(handleF)
%             end
%             
%             display(sprintf('figure saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',fName));
            
            
        end % for units
    end % for microchans
    
    
    
    
    %             %% saving stats.
    %             if exist(['./' patientID],'dir')
    %                 try
    %                     fName = sprintf('./%s/Data/%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
    %                     save([fName '.mat'],'neuronStats')
    %                 catch
    %                     mkdir(sprintf('./%s/Data/',patientID))
    %                     fName = sprintf('./%s/Data/%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
    %                     save([fName '.mat'],'neuronStats')
    %                 end
    %             elseif exist('./Data','dir')
    %                 fName = sprintf('./Data/%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
    %                 save([fName '.mat'],'neuronStats')
    %             else
    %                 fName = sprintf('%s_session_%d_CoherenceStats_alignedOn_%s',patientID,patientID,sessionNum,alignName);
    %                 save([fName '.mat'],'neuronStats')
    %             end
    
    
end % looping over channels for each align spot.

