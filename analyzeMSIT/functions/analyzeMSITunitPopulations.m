function [outputVariable] = analyzeMSITunitPopulations(patientID,sessionNum,nevFile)
%ANALYZEMSITUNITPOPULATIONS analyze MSIT unit populations
%
%   [neuronStats,Rates] = analyzeMSITunits(patientID,sessionNumber,nevFile)
%

% Author: EHS 20160712
% VersionControl: https://github.com/elliothsmith/MSIT-analysis



%% //conflict colors//
col0 = [183 30 103]./255;
col1a = [246 139 31]./255;
col1b = [0 166 81]./255;
col2 = [82 79 161]./255;


%% loading data from NEV file
display('loading data...')
% [pathstr, name, ext] = fileparts(nevFile);

ext = nevFile(end-3:end);
if strcmp(ext,'.nev')
    NEV = openNEV(nevFile,'read');
elseif strcmp(ext,'.mat')
    load(nevFile);
end


%% organizing behavioral markers
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% [20160824] removing practice trial triggers and updating number of trials.
if ~(strcmp(patientID,'CUCX2') && isequal(sessionNum,3))
    [trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
    nTrials = sum(trigs==90);
end


%% creating neural timing variable
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];
ChanUnitTimestamp(ChanUnitTimestamp(:,2)==255,:) = [];

inclChans = unique(ChanUnitTimestamp(:,1));
nChans = length(inclChans);


%% reaction time calculation
cueTimes = trigTimes(trigs>=1 & trigs<28);
responseTimes = double(trigTimes(trigs>=100 & trigs<=104));
if length(cueTimes)<length(responseTimes)
    RTs = responseTimes(1:length(cueTimes))-cueTimes;
elseif length(cueTimes)>length(responseTimes)
    RTs = responseTimes-cueTimes(1:length(responseTimes));
else
    RTs = responseTimes-cueTimes;
end


%% building Tensor:: Starting with just making a matrix over units, then separate into trials.
for ch = 1:nChans
    ChanUnits(ch) = max(unique(ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch),2)));
    % creating a little array if channel and unit indices for each "Neuron"
    if ch == 1
        unitIdx = [repmat(inclChans(ch),ChanUnits(ch),1) [1:ChanUnits(ch)]'];
    elseif ch > 1
        unitIdx = [unitIdx; repmat(inclChans(ch),ChanUnits(ch),1) [1:ChanUnits(ch)]'];
    end
end
numUnits = sum(ChanUnits);
% sanity check to make sure the previous step didn't go ostensibly wrong.
if size(unitIdx,1)~=numUnits
    display('oops. number of units does not match the index variable')
end

% How much Time is represented in the timestamps?
timeSize = ceil(max(ChanUnitTimestamp(:,3))*1000);

% builidng a binary matrix for APs.
APmat = zeros(numUnits,timeSize);

% populating that matrix
for ds = 1:size(ChanUnitTimestamp,1)
    APmat(sum(repmat(ChanUnitTimestamp(ds,1:2),numUnits,1)==unitIdx,2)==2,round(ChanUnitTimestamp(ds,3)*1000)) = 1;
end

% sanity check to make sure the previous step didn't go ostensibly wrong.
if sum(sum(APmat))~=size(ChanUnitTimestamp,1)
    display('oops. number of units does not match the index variable')
end
% removing channels without APs and updating number of units
badIdx = sum(APmat,2)==0;
APmat(badIdx,:) = [];
unitIdx(badIdx,:) = [];
%sanity Check
display(['removed unit ' num2str(find(badIdx)) ' due to low firing rate.'])
numUnits = size(APmat,1);


%% convolve APmat with a gaussian kernel to estimate firing rate.
% 1) constructing the appropriate kernel
SIGMA = 512; % standard deviation of Kernel in Milliseconds.
HSIZE = [1 1000]; % size of the gaussian kernel (1 second)
Gauss1 = fspecial('gaussian',HSIZE,SIGMA);
GaussN = repmat(Gauss1./max(Gauss1),numUnits,1);

% 2) convolving with gaussian
APrate = conv2(APmat,Gauss1,'same');


%% parsing behavior
trialType = zeros(1,nTrials);
condition = trigs(trigs>=1 & trigs<=27);


%% setting up codes for PSTHs over conflict types.
% These are the correct codes. Double Checked on 20160216
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)


%% [20160927] aligning trials into tensors for each cue and response.
for aS = 1
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            responseTimes = trigTimes(trigs>=100 & trigs<=105);
            %% timing (seconds)
            pre = 0.5;
            post = 5;
        case 2
            alignName = 'Response';
            trialStarts = trigTimes(trigs>=100 & trigs<=105);
            %% timing (seconds)
            pre = 1;
            post = 1;
    end
    
    % setting timing
    trialStartMS = ceil(trialStarts*1000) - (pre*1000);
    trialEndMS = trialStartMS + ((pre+post)*1000);
    
    % initializing varaibles.
    rateTensor = zeros(nTrials,size(APrate,1),((pre+post)*1000)+1);
    APtensor = zeros(nTrials,size(APmat,1),((pre+post)*1000)+1);
    
    %% TODO:: remove discharge trials.
    
    
    % looping over trials.
    for tt = 1:nTrials
        % stacking AP details into a tensor.
        if trialEndMS(tt)<length(APrate)
            rateTensor(tt,:,:) = APrate(:,trialStartMS(tt):trialEndMS(tt));
            APtensor(tt,:,:) = APmat(:,trialStartMS(tt):trialEndMS(tt));
        else
            leftoverData = zeros(size(APrate,1),trialEndMS(tt)-length(APrate));
            rateTensor(tt,:,:) = [APrate(:,trialStartMS(tt):end) leftoverData];
            APtensor(tt,:,:) = [APmat(:,trialStartMS(tt):end) leftoverData];
        end
        
        %% [20161011] now to examine several measures of single-trial representations across the population
        %   1) coherence across neurons.
        %   2) mutual information for AP timing.
        
        
        %% [20161011] setting up spatial coherence. pairwise?
        
        
        %% [20161011] pairwise mutual information for spike timing.
        % defining an encoding window for specific trial.
        encodingWindow{tt} = ceil(trialStarts(tt)*1000):ceil(responseTimes(tt)*1000);
        
        %         % data and time saving.
        %         dataPath = sprintf('/home/elliot/data/msit_units/%s/Data/',patientID);
        %         mutualInfoFile = sprintf('%s_pairwiseMutualInformation_MSIT_session%d.mat',patientID,sessionNum);
        %
        %         if ~exist([dataPath mutualInfoFile],'file')
        % permutations over channels - ||TAKES A LONG TIME TO RUN}||
        C = combnk(1:size(APmat,1),2);
        % looping over channel combinations
        tic
        coincidenceMeasure = 'R'; % can be 'MI' or 'CI' or 'R'
        for cmb = 1:size(C,1)
            switch coincidenceMeasure
                case {'MI'}
                    % updating user
                    if isequal(mod(cmb,1000),0)
                        display([patientID ' trial:' num2str(tt) ' calculated mutual information for combinations ' num2str(cmb) ' of ' num2str(size(C,1))])
                    end
                    % calculating mutual information between a pair of spike trains
                    Mx = APmat(C(cmb,1),encodingWindow{tt});
                    My = APmat(C(cmb,2),encodingWindow{tt});
                    % mutual info calculation
                    MI(tt,cmb) = MutualInfo(Mx,My,0);
                    %                     MInorm(tt,cmb) =
                case {'CI'}
                    % updating user
                    if isequal(mod(cmb,1000),0)
                        display([patientID ' calculated coincidence index for combinations ' num2str(cmb) ' of ' num2str(size(C,1))])
                    end
                    % calculating coincidence index
                    numCooccurringSpikes = sum(APmat(C(cmb,1),encodingWindow{tt})+APmat(C(cmb,2),encodingWindow{tt})==2);
                    numPossiblyCooccurringSpikes = max([sum(APmat(C(cmb,1),encodingWindow{tt})) sum(APmat(C(cmb,2),encodingWindow{tt}))]);
                    CI(tt,cmb) = numCooccurringSpikes/numPossiblyCooccurringSpikes;
                    
            end
        end
        A = toc;
        display(sprintf('pairwise mutual information calculation for trial %d took %d seconds.',tt,A))
        
        %             save([dataPath mutualInfoFile]);
        %         elseif isequal(exist([dataPath mutualInfoFile],'file'),2)
        %             load([dataPath mutualInfoFile])
        %         end
        switch coincidenceMeasure
            case {'R'}
                %                 timeWin = [0 median(RTs)];
                %                     sigmas = [1 2 4 8 16 32 64 128 256 512 1024 2048]; % in ms
                %                     for sw = 1:length(sigmas)
                %
                %                         sprintf('calculating spike-timing reliability for sigma = %d',sigmas(sw))
                %                         kernelWidth = sigmas(sw) ./1000;
                %                         % in case there aren't any spikes...
                %                         for tt = 1:nTrials
                %                             %% calculating single-trial firing rate.
                %                             try
                %                                 [Rtft(tt,:),tft,~] = psth(data.unit(cl).trial(tt), kernelWidth, 'n', [0 pre+post]);
                %                             catch
                %                                 Rtft(tt,:) = nan(1,length(R));
                %                             end
                %                         end
                %                         tsecft = tft-repmat(pre,1,length(tft));
                
                % calculating mean pearson correlation for each trial condition.
                %                         reliability(sw).signaWidth = sigmas(sw);
                
                R{tt} = corr(APrate(:,encodingWindow{tt}));
                
        end
    end
    
    quantVis = 'linreg'; % can be either 'bar' or 'linreg'
    switch coincidenceMeasure
        case {'MI'}
            % Taking the median value from across the distribution
            MIbar = nanmedian(MI,2);
            
            switch quantVis
                case {'bar'}
                    % plotting mean MI results
                    figure
                    hold on
                    bar(1:4,[mean(MIbar(trialType==1)) mean(MIbar(trialType==2)) mean(MIbar(trialType==3)) mean(MIbar(trialType==4))])
                    errorbar(1:4,[mean(MIbar(trialType==1)) mean(MIbar(trialType==2)) mean(MIbar(trialType==3)) mean(MIbar(trialType==4))],...
                        [std(MIbar(trialType==1))./sum(trialType==1) std(MIbar(trialType==2))./sum(trialType==2) std(MIbar(trialType==3))./sum(trialType==3) std(MIbar(trialType==4))./sum(trialType==4)],...
                        'marker','none','linestyle','none','linewidth',3,'color','k')
                    hold off
                    axis square
                    ylabel('median pairwise mutual information (bits)','fontsize',14)
                    set(gca,'fontsize',14,'linewidth',1.5,'XtickLabel',{'','None','Simon','Eriksen','Both',''})
                    
                    outputVariable = MI;
                case {'linreg'}
                    % linear regression
                    %                     keyboard
                    [rhoMI,pMI] = corr(RTs',MIbar,'type','spearman','tail','both');
                    
                    % making a color map
                    cMap = zeros(length(RTs),3);
                    cMap(trialType==1,:) = repmat(col0,sum(trialType==1),1);
                    cMap(trialType==2,:) = repmat(col1a,sum(trialType==2),1);
                    cMap(trialType==3,:) = repmat(col1b,sum(trialType==3),1);
                    cMap(trialType==4,:) = repmat(col2,sum(trialType==4),1);
                    %% ooh: could make corcle sizes another variable:: theta coh?
                    
                    
                    % scatter plot of RTs and MI over trials
                    figure(112358)
                    hold on
                    scatter(RTs',MIbar,repmat(30,length(RTs),1),cMap,'filled')
                    text(min(RTs),max(MIbar),['rho = ' num2str(rhoMI) ', p = ' num2str(pMI)])
                    hold off
                    set(gca, 'linewidth',2,'fontsize',14)
                    
                    %                     plotOverSessions = 1;
                    %                     if ~plotIOverSessions
                    %                         saveas()
                    %                         close
                    %                     else
                    %                         saveas()
                    %                     end
            end
            
            outputVariable = MI;
            
        case {'CI'}
            
            % Taking the median value from across the distribution
            MIbar = nanmean(CI,2);
            
            switch quantVis
                case {'bar'}
                    % plotting mean MI results
                    figure
                    hold on
                    bar(1:4,[mean(MIbar(trialType==1)) mean(MIbar(trialType==2)) mean(MIbar(trialType==3)) mean(MIbar(trialType==4))])
                    errorbar(1:4,[mean(MIbar(trialType==1)) mean(MIbar(trialType==2)) mean(MIbar(trialType==3)) mean(MIbar(trialType==4))],...
                        [std(MIbar(trialType==1))./sum(trialType==1) std(MIbar(trialType==2))./sum(trialType==2) std(MIbar(trialType==3))./sum(trialType==3) std(MIbar(trialType==4))./sum(trialType==4)],...
                        'marker','none','linestyle','none','linewidth',3,'color','k')
                    hold off
                    axis square
                    ylabel('median pairwise mutual information (bits)','fontsize',14)
                    set(gca,'fontsize',14,'linewidth',1.5,'XtickLabel',{'','None','Simon','Eriksen','Both',''})
                    
                case {'linreg'}
                    % linear regression
                    %                     keyboard
                    [rhoMI,pMI] = corr(RTs',MIbar,'type','spearman','tail','both');
                    
                    % making a color map
                    cMap = zeros(length(RTs),3);
                    cMap(trialType==1,:) = repmat(col0,sum(trialType==1),1);
                    cMap(trialType==2,:) = repmat(col1a,sum(trialType==2),1);
                    cMap(trialType==3,:) = repmat(col1b,sum(trialType==3),1);
                    cMap(trialType==4,:) = repmat(col2,sum(trialType==4),1);
                    %% ooh: could make corcle sizes another variable:: theta coh?
                    
                    
                    % scatter plot of RTs and MI over trials
                    figure(112358)
                    hold on
                    scatter(RTs',MIbar,repmat(30,length(RTs),1),cMap,'filled')
                    text(min(RTs),max(MIbar),['rho = ' num2str(rhoMI) ', p = ' num2str(pMI)])
                    hold off
                    set(gca, 'linewidth',2,'fontsize',14)
                    
                    plotOverSessions = 0;
                    if plotOverSessions
                        saveas(112358,'/home/elliot/data/msit_units/CUCX2/Figs/populations/spikeTimingCoincidence_all3sessions.pdf')
                    else
                        saveas(112358,['/home/elliot/data/msit_units/CUCX2/Figs/populations/spikeTimingCoincidence_session' num2str(sessionNum) '.pdf'])
                        close (112358)
                    end
            end
            
            outputVariable = CI;
        case {'R'}
            switch quantVis
                case {'linreg'}
                    % making a color map
                    cMap = zeros(length(RTs),3);
                    cMap(trialType==1,:) = repmat(col0,sum(trialType==1),1);
                    cMap(trialType==2,:) = repmat(col1a,sum(trialType==2),1);
                    cMap(trialType==3,:) = repmat(col1b,sum(trialType==3),1);
                    cMap(trialType==4,:) = repmat(col2,sum(trialType==4),1);
                    
                    % scatter plot of reliability vs RT
                    figure(112358)
                    hold on
                    scatter(RTs',R,repmat(30,length(RTs),1),cMap,'filled')
                    hold off
                    set(gca, 'linewidth',2,'fontsize',14)
                    
                    
                    saveas(112358,['/home/elliot/data/msit_units/CUCX2/Figs/populations/spikeTimingReliability_perTrial_session' num2str(sessionNum) '.pdf'])
                    close (112358)
            end
    end
    
    
    %     MIbar = nanmean(zscore(MI(:,1:end),[],2));
    %     tsecMI = tsecSpikes(1:length(MIbar));
    
end









%
%
%
% %% Rasters and PSTHs
% for aS = 1:length(data)
%     % which alignment spot
%     switch aS
%         case 1
%             alignName = 'Cue';
%             trialStarts =  trigTimes(trigs>=1 & trigs<28);
%             %% timing (seconds)
%             pre = 2;
%             post = 3;
%         case 2
%             alignName = 'Response';
%             trialStarts =  trigTimes(trigs>=100 & trigs<=105);
%             %% timing (seconds)
%             pre = 3;
%             post = 3;
%     end
%
%     % looping over units
%     for ch = 1:size(data(aS).channel,2)
%         for un = 1:size(data(aS).channel(ch).unit,2)
%             for tt = 1:nTrials
%                 %% PSTHs
%                 % calculating psths
%                 kernelWidth = 25  ./1000;
%                 [R(tt,:),t,E] = psth(data(aS).channel(ch).unit(un).trial(), kernelWidth, 'n', [0 pre+post]);
%                 tsec = t-repmat(pre,1,length(t));
%
%             end
%             keyboard
%
%         end % looping over units for each channel and align spot
%     end % looping over channels for each align spot.
%
%
%     %% saving figures.
%     figFlag = 1;
%     if figFlag
%         if exist('./Figs','dir')
%             if exist('./Figs/FiringRate/','dir')
%                 fName = sprintf('./Figs/FiringRate/%s_session_%d_Channel_%d_Unit_%d_Conflict_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
%                 saveas(aS,fName, 'pdf')
%                 close(aS)
%             else
%                 mkdir('./Figs/FiringRate/')
%                 fName = sprintf('./Figs/FiringRate/%s_session_%d_Channel_%d_Unit_%d_Conflict_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
%                 saveas(aS,fName, 'pdf')
%                 close(aS)
%             end
%         else
%             try
%                 fName = sprintf('./Figs/%s_session_%d_Channel_%d_Unit_%d_Conflict_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
%                 saveas(aS,fName, 'pdf')
%                 close(aS)
%             catch
%                 fName = sprintf('%s_session_%d_Channel_%d_Unit_%d_Conflict_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
%                 saveas(aS,fName, 'pdf')
%                 close(aS)
%             end
%         end
%     end
%
%
%     %% saving stats.
%     if exist('./Data','dir')
%         fName = sprintf('./Data/%s_session_%d_ConflictStats_alignedOn_%s',patientID,sessionNum,alignName);
%         save([fName '.mat'],'neuronStats')
%     elseif ~exist('./Data','dir')
%         mkdir('./Data/')
%         fName = sprintf('./Data/%s_session_%d_ConflictStats_alignedOn_%s',patientID,sessionNum,alignName);
%         save([fName '.mat'],'neuronStats')
%     else
%         fName = sprintf('%s_session_%d_ConflictStats_alignedOn_%s',patientID,sessionNum,alignName);
%         save([fName '.mat'],'neuronStats')
%     end
%
% end % looping over align spots (Stimulus & response)
%
