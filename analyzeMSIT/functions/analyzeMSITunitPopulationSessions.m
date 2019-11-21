function [outputVariable] = analyzeMSITunitPopulationSessions(patientID,startSession)
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


%% TODO: get patient Directory.
if isequal(exist(['./' patientID],'dir'),7)
    display('patient directory exists in cwd.')
    patientDir = './';
else
    patientDir = uigetdir('~/','Select directory with patient sub-directories');
    cd(patientDir)
end


Rsessions = [];
RTsessions = [];
trialTypeSessions = [];
%% running the functions from the input arguments.
if ~isequal(exist(patientID,'dir'),7)
    display('no patient directory found')
else
    % entering patient directory`
    cd(fullfile(patientDir,patientID))
    fullPath = pwd;
    
    % creating data directory
    if ~exist('Data','dir')
        display(sprintf('making data directory. You should probably put your data in %s',strcat(fullPath,'Data')))
        mkdir(fullPath,'Data');
        if isunix
            addpath('./Data')
        else
            addpath('.\Data')
        end
    end
    % creating figure directory
    if ~exist('Figs','dir')
        mkdir(fullPath,'Figs');
    end
    
    % finding .nev files
    dirList = dir(fullfile(fullPath,'Data','*.nev'));
    display(['found ' num2str(length(dirList)) ' NEV files!'])
    if isequal(length(dirList),0)
        display(sprintf('\t... now looking for .mattified NEVs... \nFormatting note: this code will look for files named *NEV.mat. \nIf you have .mattified NEV files without "NEV" at the end of their name, consider renaming.'))
        dirList = dir(fullfile(fullPath,'Data','*NEV.mat'));
        display(['found ' num2str(length(dirList)) ' .mattified NEV files!'])
    end
    
    
    
    for sessionNum = 1+startSession:length(dirList)
        nevFile = fullfile(fullPath,'Data',dirList(sessionNum).name);
        
        %% loading data from NEV file
        display(['loading data for session ' num2str(sessionNum)])
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
        
        clear APrate
        %% looping over kernel widths
        sigmas = [1 2 4 8 16 32 64 128 256 512 1024 2048]; % in ms
        for sw = 1:length(sigmas)
            % updating user
            display(sprintf('calculating firing rate across the full recording with for %d s.d. kernel',sigmas(sw)))
            
            %% convolve APmat with a gaussian kernel to estimate firing rate.
            % 1) constructing the appropriate kernel
            SIGMA = sigmas(sw); % standard deviation of Kernel in Milliseconds.
            HSIZE = [1 1000]; % size of the gaussian kernel (1 second)
            Gauss1 = fspecial('gaussian',HSIZE,SIGMA);
            GaussN = repmat(Gauss1./max(Gauss1),numUnits,1);
            
            % 2) convolving with gaussian
            APrate(:,:,sw) = conv2(APmat,Gauss1,'same');
        end
        
        %% parsing behavior
        badBehavior = (RTs<0.3 | RTs>5);
        trialType = zeros(1,nTrials);
        condition = trigs(trigs>=1 & trigs<=27);
        
        
        %% setting up codes for PSTHs over conflict types.
        % These are the correct codes. Double Checked on 20160216
        trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
        trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
        trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
        trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)
        
        
        %% removing bad trials
        %         keyboard
        RTs (badBehavior) = [];
        trialType (badBehavior) = [];
        
        RTcell{sessionNum} = RTs';
        trialTypeCell{sessionNum} = trialType;
        
        nTrials = length(RTs);
        
        
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
                display(['trial: ' num2str(tt)])
                %                 % stacking AP details into a tensor.
                %                 if trialEndMS(tt)<length(APrate)
                %                     rateTensor(tt,:,:) = APrate(:,trialStartMS(tt):trialEndMS(tt));
                %                     APtensor(tt,:,:) = APmat(:,trialStartMS(tt):trialEndMS(tt));
                %                 else
                %                     leftoverData = zeros(size(APrate,1),trialEndMS(tt)-length(APrate));
                %                     rateTensor(tt,:,:) = [APrate(:,trialStartMS(tt):end) leftoverData];
                %                     APtensor(tt,:,:) = [APmat(:,trialStartMS(tt):end) leftoverData];
                %                 end
                
                %% [20161011] now to examine several measures of single-trial representations across the population
                %   1) coherence across neurons.
                %   2) mutual information for AP timing.
                
                
                %% [20161011] setting up spatial coherence. pairwise?
                
                
                %% [20161011] pairwise mutual information for spike timing.
                encWin = '4secs'; % can be 'meanRT' or 'eachRT' or '4secs'
                if strcmp(encWin,'eachRT')
                    % defining an encoding window for specific trial.
                    encodingWindow{tt} = ceil(trialStarts(tt)*1000):ceil(responseTimes(tt)*1000);
                elseif strcmp(encWin,'meanRT')
                    % defining a constant encoing window == mean RT
                    encodingWindow{tt} = ceil(trialStarts(tt)*1000):ceil((trialStarts(tt)+median(RTs))*1000);
                elseif strcmp(encWin,'4secs')
                    encodingWindow{tt} = ceil(trialStarts(tt)*1000:(trialStarts(tt)+5)*1000);
                end
                
                
                %% picking the measure of spike timing coincidence
                coincidenceMeasure = 'MI'; % can be 'MI' or 'CI' or 'R'
                normalizationOperation = 'mean';
                
                %% data and time saving
                dataPath = sprintf('/home/elliot/data/msit_units/CUCX2/Data/populations/');
                mutualInfoFile = sprintf('%s_%sdata_%s_%s.mat',patientID,coincidenceMeasure,normalizationOperation,encWin);
                if ~exist([dataPath mutualInfoFile],'file')
                    
                    % permutations over channels - ||TAKES A LONG TIME TO RUN}||
                    C = combnk(1:size(APmat,1),2);
                    
                    
                    %% calculating mutual information over combinations. 20161028::data should be loaded
                    % looping over channel combinations
                    for cmb = 1:size(C,1)
                        if isequal(coincidenceMeasure,'MI')
                            
                            % updating user
                            if isequal(mod(cmb,1000),0)
                                display([patientID ' trial:' num2str(tt) ' calculated mutual information for combinations ' num2str(cmb) ' of ' num2str(size(C,1))])
                            end
                            % calculating mutual information between a pair of spike trains
                            Mx = APmat(C(cmb,1),encodingWindow{tt});
                            My = APmat(C(cmb,2),encodingWindow{tt});
                            % mutual info calculation
                            if (isequal(sum(Mx),0))
                                MI{sessionNum}(tt,cmb) = 0;
                            else
                                MI{sessionNum}(tt,cmb) = MutualInfo(Mx,My,0);
                            end
                        end
                    end
                    
                    
                    %% calculating coincidence index
                    if isequal(coincidenceMeasure,'CI')
                        % updating user
                        disp(' calculated coincidence index for ...')
                        %                         keyboard
                        % calculating coincidence index
                        meanCooccurringSpikes(tt) = nanmean(nansum(APmat(:,encodingWindow{tt})));
                        totalSpikes(tt) = nansum(nansum(APmat(:,encodingWindow{tt})));
                    elseif strcmp(coincidenceMeasure,'R')
                        for sw2 = 1:length(sigmas)
                            Rmat = corr(APrate(:,encodingWindow{tt},sw2));
                            R(tt,sw2) = nanmean(Rmat(:));
                            clear Rmat
                        end
                    end
                end
            end
            %             keyboard
            %             CI{sessionNum} = (meanCooccurringSpikes./totalSpikes);
            
            %
            %             A = toc;
            %             display(sprintf('pairwise mutual information calculation took %d seconds for all trials.',A))
            
            
        end % cue aligned
        
        if strcmp(coincidenceMeasure,'R')
            Rsessions = cat(1,Rsessions,R);
            clear R
            RTsessions = cat(2,RTsessions,RTs);
            trialTypeSessions = cat(2,trialType,trialTypeSessions);
        end
        
    end % looping over sessions
    quantVis = 'bar'; % can be either 'bar' or 'linreg'
    
    
    %% now plotting and quantifying
    switch coincidenceMeasure
        case {'MI'}
            
            % saving or loading data if MI was calculated this time around.
            if ~exist([dataPath mutualInfoFile],'file')
                % measure of centrality over the population
                MI1 = nanmedian(MI{1},2);
                MI2 = nanmedian(MI{2},2);
                MI3 = nanmedian(MI{3},2);
                
                % calculating outlier indices
                oMI1 = true(size(MI1));
                oMI2 = true(size(MI2));
                oMI3 = true(size(MI3));
                
                % removing outlier indices
                oMI1(outliers(MI1)) = false;
                oMI2(outliers(MI2)) = false;
                oMI3(outliers(MI3)) = false;
                
                % Taking the median value from across the distribution
                MIbar = cat(1,nanzscore(MI1(oMI1)),nanzscore(MI2(oMI2)),nanzscore(MI3(oMI3)));
                allRTs = cat(1,RTcell{1}(oMI1),RTcell{2}(oMI2),RTcell{3}(oMI3));
                allTrialTypes = cat(1,trialTypeCell{1}(oMI1)',trialTypeCell{2}(oMI2)',trialTypeCell{3}(oMI3)');
                
                save([dataPath mutualInfoFile],'MIbar','allRTs','allTrialTypes','MI','RTcell','trialTypeCell')
                
            else
                load([dataPath mutualInfoFile])
            end
            
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
                    keyboard
                    [rhoMI,pMI] = corr(allRTs,MIbar,'type','spearman','tail','both');
                    
                    % making a color map
                    cMap = zeros(length(allRTs),3);
                    cMap(allTrialTypes==1,:) = repmat(col0,sum(allTrialTypes==1),1);
                    cMap(allTrialTypes==2,:) = repmat(col1a,sum(allTrialTypes==2),1);
                    cMap(allTrialTypes==3,:) = repmat(col1b,sum(allTrialTypes==3),1);
                    cMap(allTrialTypes==4,:) = repmat(col2,sum(allTrialTypes==4),1);
                    %% ooh: could make circle sizes another variable:: theta coh?
                    
                    
                    % scatter plot of RTs and MI over trials
                    figure(112358)
                    hold on
                    scatter(allRTs,MIbar,repmat(30,length(allRTs),1),cMap,'filled')
                    text(min(allRTs),max(MIbar),['rho = ' num2str(rhoMI) ', p = ' num2str(pMI)])
                    hold off
                    set(gca, 'linewidth',2,'fontsize',14)
                    
                    saveas(112358,'/home/elliot/data/msit_units/CUCX2/Figs/populations/spikeTimingCoincidence_mutualInformation_all3sessions.pdf')
                    saveas(112358,'/home/elliot/Dropbox/Figs/MSIT-ACC-dlPFC/spikeTimingCoincidence_mutualInformation_all3sessions.pdf')
            end
            
            outputVariable = MI;
            
        case {'CI'}
            
            MIbar = zscore(cat(1,CI{1}',CI{2}',CI{3}'));
            allRTs = cat(1,RTcell{1},RTcell{2},RTcell{3});
            allTrialTypes = cat(1,trialTypeCell{1}',trialTypeCell{2}',trialTypeCell{3}');
            
            
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
                    keyboard
                    [rhoMI,pMI] = corr(allRTs,MIbar,'type','spearman','tail','both');
                    
                    % making a color map
                    cMap = zeros(length(allRTs),3);
                    cMap(allTrialTypes==1,:) = repmat(col0,sum(allTrialTypes==1),1);
                    cMap(allTrialTypes==2,:) = repmat(col1a,sum(allTrialTypes==2),1);
                    cMap(allTrialTypes==3,:) = repmat(col1b,sum(allTrialTypes==3),1);
                    cMap(allTrialTypes==4,:) = repmat(col2,sum(allTrialTypes==4),1);
                    %% ooh: could make circle sizes another variable:: theta coh?
                    
                    %                     keyboard
                    % scatter plot of RTs and MI over trials
                    figure(112358)
                    hold on
                    scatter(allRTs,MIbar,repmat(30,length(allRTs),1),cMap,'filled')
                    text(min(allRTs),max(MIbar),['rho = ' num2str(rhoMI) ', p = ' num2str(pMI)])
                    hold off
                    set(gca, 'linewidth',2,'fontsize',14)
                    
                    
                    saveas(112358,'/home/elliot/Dropbox/Figs/MSIT-ACC-dlPFC/spikeTimingCoincidence_all3sessions.pdf')
                    saveas(112358,'/home/elliot/data/msit_units/CUCX2/Figs/populations/spikeTimingCoincidence_all3sessions.pdf')
            end
            
            outputVariable = CI;
            
        case {'R'}
            switch quantVis
                case {'linreg'}
                    % making a color map fro trial types
                    cMap = zeros(length(trialTypeSessions),3);
                    cMap(trialTypeSessions==1,:) = repmat(col0,sum(trialTypeSessions==1),1);
                    cMap(trialTypeSessions==2,:) = repmat(col1a,sum(trialTypeSessions==2),1);
                    cMap(trialTypeSessions==3,:) = repmat(col1b,sum(trialTypeSessions==3),1);
                    cMap(trialTypeSessions==4,:) = repmat(col2,sum(trialTypeSessions==4),1);
                    
                    
                    %                     load('spikeTimingReliabilityData_allKernels')
                    
                    keyboard
                    for sw3 = 1:length(sigmas)
                        % scatter plot of reliability vs RT
                        figure(sw3)
                        hold on
                        % now doing and plotting linear regression
                        %                     regWin
                        mdl = fitlm(RTsessions,Rsessions(:,sw3),'linear') %,'Intercept',true,'RobustOpts','on'
                        plot(linspace(-1,3,10),repmat(mdl.Coefficients{2,1},1,10).*linspace(-1,3,10)+repmat(mdl.Coefficients{1,1},1,10),'--k','linewidth',2);
                        
                        % plotting data
                        scatter(RTsessions,Rsessions(:,sw3),repmat(10,length(RTsessions),1),cMap,'filled')
                        hold off
                        set(gca, 'linewidth',2,'fontsize',14)
                        ylim ([0 1])
                        
                        saveas(sw3,['/home/elliot/Dropbox/Figs/MSIT_ACC_dlPFC/spikeTimingReliability_perTrial_allsessions_kernelSD' num2str(sigmas(sw3)) '.pdf'])
                        %                                              close (112358)
                    end
            end
            
    end
    
end % no patient directory catch.






