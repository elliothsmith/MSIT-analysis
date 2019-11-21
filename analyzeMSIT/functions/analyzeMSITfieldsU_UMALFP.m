function [fieldStats] = analyzeMSITfieldsU_UMALFP(patientID,sessionNum,nevFile,tf)
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
nsFile = fullfile(dataPath,[nvName(1:19) '.ns5']);
if ~exist(nsFile,'file')
    nsList = dir(fullfile(dataPath,[nvName(1:8) '*.ns5']));
    nsFile = fullfile(dataPath,nsList(1).name);
end


% parsing files and loading data.
if strcmp(nvExt,'.nev')
    NEV = openNEV(nevFile,'read');
%     if isequal(exist(nsFile,'file'),0)
%         [nsFile,nsPath,~] = uigetfile('*.ns*','Select the correct NSx file');
%         NS3 = openNSx(fullfile(nsPath,nsFile));
%     else
%         NS3 = openNSx(nsFile);
%     end
elseif strcmp(nvExt,'.mat')
    load(nevFile);
%     if isequal(exist(nsFile,'file'),0)
%         [nsFile,nsPath,~] = uigetfile('*.ns*','Select the correct NSx file');
%         NS3 = openNSx(fullfile(nsPath,nsFile));
%     else
%         NS3 = openNSx(nsFile);
%     end
end


%% loading LFP
load([nsFile(1:end-4) '_downsampledMeanUMALFP.mat'])


%% organizing important task parameters.
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% [20160824] removing practice trial triggers and updating number of trials.
[trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
nTrials = sum(trigs==90);


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
responses = double(trigs(trigs>=100 & trigs<=104))';
RTs = responses(1:length(cues))-cues;



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
nTrials = 69;
trialType(nTrials+1:end) = [];


% %% [20160610] denoising LFP
% denoiseMethod = 'notch'
% switch denoiseMethod
%     case 'PCA'
%         display('denoising using PCA...')
%         dData = remove1stPC(double(NS3.Data));
%         display('...done.')
%     case 'notch'
%         Wo = 60/(params.Fs/2);  BW = Wo/50;
%         [b,a] = iirnotch(Wo,BW);
%         %        freqz(b,a);
%         for c = 1:size(NS3.Data{2},1)
%             display(sprintf('applying notch filter to channel %d',c))
%             dData(c,:) = filtfilt(b,a,double(NS3.Data{2}(c,:)));
%         end
% end


%% building LFP tensor and plotting ERPs and spectrograms.
display('Aligning LFP data on stimulus and response.');
% Channel labels
LFPlabels = 'UMALFP'; % anterior ECoG contact to the array.
% cue timing
pre = 2;
post = 3;
% initializing LFPmat
LFPmat = zeros(((pre+post)*params.Fs)+1,length(goods),length(LFPlabels),2); % 2 for cue and response aligned.
% for loop to save multiple epochs
for aS = 1
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            %             nTrials = sum(trigs>=1 & trigs<28);
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>=100 & trigs<=105);
            %             trigTimes(trigs>=100 & trigs<104)
            % response timing
            pre = 3;
            post = 2;
    end
    
    %% time space for LFP samples.
    tsec = linspace(-pre,post,size(LFPmat,1));
    
    % looping over channels.
    for ch = 1:length(LFPlabels)
        % looping over trials.
        for tt = 1:nTrials
            
            %% constructing LFP tensor.
            try
                LFPmat(:,tt,ch,aS) = dsData(floor((trialStarts(goods(tt))-pre)*params.Fs):floor((trialStarts(goods(tt))+post)*params.Fs));
            catch
                LFPmat(:,tt,ch,aS) = dsData(ceil((trialStarts(goods(tt))-pre)*params.Fs):floor((trialStarts(goods(tt))+post)*params.Fs));
                
            end
            
            %% [20160610] calculating single trial spectrograms here IF THEY ARE NEEDED.
            params.trialave = 0;
            if isequal(mod(tt,50),0)
                display(sprintf('calculated spectrograms for for first %d trials out of %d...',tt,nTrials))
            end
            [Stf(:,:,tt),t,f] = mtspecgramc(LFPmat(:,tt,ch,aS),movingWin,params);
            tspec = linspace(t(1)-pre,t(end)-pre,length(t));
            
            % averaging to collect spectra.
            Sf_encodingWindow(:,tt) = mean(Stf(tspec>0 & tspec<RTs(tt),:,tt));
            Sf_baseline(:,tt) = mean(Stf(tspec>-1.5 & tspec<-0.5,:,tt));
            
            %             %% calculating spectra for each trial
            %             if isequal(tf,0)
            %                 params.pad = 6; %
            %                 if mod(tt,50)==0
            %                     display(sprintf('calculated trial averaged spectra for %d trials',tt))
            %                 end
            %                 [Sf(:,tt),f] = mtspectrumc(LFPmat(tsec>0 & tsec<RTs(tt),tt,ch,aS),params);
            %             end
        end % looping over trials
        
        %% calculating and plotting trial averaged spectrograms.
        switch tf
            case 1
                spectralType = 'Spectrograms';
                
                % calculating spectrograms
                params.trialave = 1;
                display('calculating trial averaged spectrogram for no conlfict trials...')
                [Stf0,~,~] = mtspecgramc(LFPmat(:,trialType==1,ch,aS),movingWin,params);
                display('calculating trial averaged spectrogram for simon conlfict trials...')
                [Stf1a,~,~] = mtspecgramc(LFPmat(:,trialType==2,ch,aS),movingWin,params);
                display('calculating trial averaged spectrogram for flanker conlfict trials...')
                [Stf1b,~,~] = mtspecgramc(LFPmat(:,trialType==3,ch,aS),movingWin,params);
                display('calculating trial averaged spectrogram for both conlfict trials...')
                [Stf2,t,f] = mtspecgramc(LFPmat(:,trialType==4,ch,aS),movingWin,params);
                
                
                %% time vectors
                tsec = linspace(-pre,post,size(LFPmat,1));
                tspec = linspace(t(1)-pre,t(end)-pre,length(t));
                
                normalizationType = 'baseline';
                switch normalizationType
                    case {'baseline'}
                        % mean spectrogram for all trials to determine baseline.
                        [Sall,t,f] = mtspecgramc(LFPmat(:,:,ch,aS),movingWin,params);
                        baselinePeriod = [-1.5 -0.5];
                        baselineSpectrum = mean(Sall(tspec>baselinePeriod(1) & tspec<baselinePeriod(2),:));
                        
                        % dividing by baseline spectrum.
                        Sn0 = transpose(Stf0./repmat(baselineSpectrum,size(Stf0,1),1));
                        Sn1a = transpose(Stf1a./repmat(baselineSpectrum,size(Stf1a,1),1));
                        Sn1b = transpose(Stf1b./repmat(baselineSpectrum,size(Stf1b,1),1));
                        Sn2 = transpose(Stf2./repmat(baselineSpectrum,size(Stf2,1),1));
                        
                    case {'frequencyband'}
                        % normalizing spectrograms by frequecny band.
                        Sn0 = normlogspec(Stf0)';
                        Sn1a = normlogspec(Stf1a)';
                        Sn1b = normlogspec(Stf1b)';
                        Sn2 = normlogspec(Stf2)';
                        
                end
                
                
                %% plotting LFP and spectrogram figures.
                handleF = ch*100;
                brdr = 0.04;
                figure(handleF)
                % first, ERPs
                plotmultipleaxes(2,1,4,brdr,handleF)
                hold on
                plot(tsec,nanmean(LFPmat(:,trialType==1,ch,aS),2),'color',col0)
                plot(tsec,nanmean(LFPmat(:,trialType==2,ch,aS),2),'color',col1a)
                plot(tsec,nanmean(LFPmat(:,trialType==3,ch,aS),2),'color',col1b)
                plot(tsec,nanmean(LFPmat(:,trialType==4,ch,aS),2),'color',col2)
                hold off
                xlim([-1 2])
                set(gca,'linewidth',2,'fontsize',16)
                title(deblank(LFPlabels))
                
                % then, spectrograms.
                plotmultipleaxes(2,4,2,brdr,handleF)
                imagesc(tspec,f,Sn0)
                xlim([-1 2])
                set(gca,'linewidth',2,'fontsize',16)
                axis xy square
                
                plotmultipleaxes(4,4,2,brdr,handleF)
                imagesc(tspec,f,Sn1a)
                xlim([-1 2])
                set(gca,'linewidth',2,'fontsize',16)
                axis xy square
                
                plotmultipleaxes(6,4,2,brdr,handleF)
                imagesc(tspec,f,Sn1b)
                xlim([-1 2])
                set(gca,'linewidth',2,'fontsize',16)
                axis xy square
                
                plotmultipleaxes(8,4,2,brdr,handleF)
                imagesc(tspec,f,Sn2)
                xlim([-1 2])
                set(gca,'linewidth',2,'fontsize',16)
                axis xy square
                
                colorbar('NorthOutside')
                colormap(jet)
                
            case 0
                spectralType = 'Spectra';
                
                Sf_baseline = nanmean(Sf_baseline,2);
                
                % collecting mean spectra
                Sf0dB = nanmean(Sf_encodingWindow(:,trialType==1)./repmat(Sf_baseline,1,sum(trialType==1)),2);
                Sf1adB = nanmean(Sf_encodingWindow(:,trialType==2)./repmat(Sf_baseline,1,sum(trialType==2)),2);
                Sf1bdB = nanmean(Sf_encodingWindow(:,trialType==3)./repmat(Sf_baseline,1,sum(trialType==3)),2);
                Sf2dB = nanmean(Sf_encodingWindow(:,trialType==4)./repmat(Sf_baseline,1,sum(trialType==4)),2);
                
                % calculating error over trials.
                Sf0err = nanstd(Sf_encodingWindow(:,trialType==1)./repmat(Sf_baseline,1,sum(trialType==1)),[],2)./sqrt(sum(trialType==1));
                Sf1aerr = nanstd(Sf_encodingWindow(:,trialType==2)./repmat(Sf_baseline,1,sum(trialType==2)),[],2)./sqrt(sum(trialType==2));
                Sf1berr = nanstd(Sf_encodingWindow(:,trialType==3)./repmat(Sf_baseline,1,sum(trialType==3)),[],2)./sqrt(sum(trialType==3));
                Sf2err = nanstd(Sf_encodingWindow(:,trialType==4)./repmat(Sf_baseline,1,sum(trialType==4)),[],2)./sqrt(sum(trialType==4));
                
                
                % setting up factors for 2-way ANOVA
                for fq = 1:length(f)
                    if f(fq)<4
                        freqs{fq} = 'delta';
                    elseif f(fq) >=4 && f(fq) <8
                        freqs{fq} = 'theta';
                    elseif f(fq) >=8 && f(fq) <12
                        freqs{fq} = 'alpha';
                    elseif f(fq) >=12 && f(fq) <25
                        freqs{fq} = 'beta';
                    elseif f(fq) >=25 && f(fq) <60
                        freqs{fq} = 'gamma';
                    elseif f(fq) >=60 && f(fq) <params.fpass(2)
                        freqs{fq} = 'high gamma';
                    end
                    % % for looking at each bin as a factor.
                    % freqs{fq} = num2str(flog(fq));
                end
                
                % reshaping the grouping variables.
                classesVec = repmat(trialType',length(freqs),1)';
                freqsVec = repmat(freqs,1,length(trialType));
                
                % reshaping S matrix
                Sfnorm_encodingWindow = Sf_encodingWindow./repmat(Sf_baseline,1,size(Sf_encodingWindow,2));
                SVec = reshape(Sfnorm_encodingWindow,1,length(freqsVec));
                
%                 if sum(isnan(SVec))~=length(SVec)
%                     
%                     
%                     %% [20160518] doing statistics on spectra.
%                     figure(ch)
%                     [p,table,stats] = anovan(SVec,{classesVec freqsVec},'varnames',{'trial types','frequency bands'},'display','off','model','interaction');
%                     %             [p,table,stats] = anovan(SVec,{freqsVec classesVec},'varnames',{'frequency bands','trial types'},'display','off','model','interaction');
%                     figure(ch*10)
%                     [comp,means,H,Gnames] = multcompare(stats,'alpha',0.01,'dimension',[1 2]);
%                     %             [Gnames  num2cell(means)]
%                     %             [Gnames([2 8 14 20],:)  num2cell(means([2 8 14 20],:))]
%                     
%                     
%                     %% saving multiple comparison Figures as .figs
%                     maximize(ch*10)
%                     pause(3)
%                     if exist('./Figs','dir')
%                         try
%                             fName = sprintf('./Figs/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_multipleComparisons',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
%                             saveas(ch*10,fName,'fig')
%                             close(ch*10)
%                         catch
%                             mkdir('./Figs/LFP/')
%                             fName = sprintf('./Figs/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_multipleComparisons',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
%                             saveas(ch*10,fName,'fig')
%                             close(ch*10)
%                         end
%                     elseif exist('./Figs/LFP/','dir')
%                         fName = sprintf('../../Figs/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_multipleComparisons',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
%                         saveas(ch*10,fName,'fig')
%                         close(ch*10)
%                     else
%                         fName = sprintf('%s_session_%d_Channel_%s_LFP%s_%sAligned_multipleComparisons',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
%                         saveas(ch*10,fName,'fig')
%                         close(ch*10)
%                     end
%                     
%                     display(sprintf('figure saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',fName));
%                     
%                     
%                     %% saving multiple comparison data
%                     fieldStats.twoWayANOVA.table = table;
%                     fieldStats.twoWayANOVA.stats = stats;
%                     fieldStats.postHocTest.comparison = comp;
%                     fieldStats.postHocTest.groupNames = Gnames
%                     fieldStats.postHocTest.meansAndError = means;
%                     
%                     if exist('./Data','dir')
%                         try
%                             fName = sprintf('./Data/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_multipleComparisons',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
%                             save([fName '.mat'],'fieldStats')
%                         catch
%                             mkdir('./Data/LFP/')
%                             fName = sprintf('./Data/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_multipleComparisons',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
%                             save([fName '.mat'],'fieldStats')
%                         end
%                     elseif exist('./Data/LFP/','dir')
%                         fName = sprintf('../../Data/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_multipleComparisons',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
%                         save([fName '.mat'],'fieldStats')
%                     else
%                         fName = sprintf('%s_session_%d_Channel_%s_LFP%s_%sAligned_multipleComparisons',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
%                         save([fName '.mat'],'fieldStats')
%                     end
%                     display(sprintf('LFP stats data saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',fName));
%                     
%                 end % in case of NaNs, skip stats...
%                 
                
                
                %% plotting spectra
                handleF = ch*100;
                brdr = 0.04;
                figure(handleF)
                % first, ERPs
                plotmultipleaxes(1,1,2,brdr,handleF)
                hold on
                plot(tsec,nanmean(LFPmat(:,trialType==1,ch,aS),2),'color',col0)
                plot(tsec,nanmean(LFPmat(:,trialType==2,ch,aS),2),'color',col1a)
                plot(tsec,nanmean(LFPmat(:,trialType==3,ch,aS),2),'color',col1b)
                plot(tsec,nanmean(LFPmat(:,trialType==4,ch,aS),2),'color',col2)
                text(0,40,['mean RT: ' num2str(mean(RTs(trialType==1))) ' s'],'color',col0,'fontsize',16,'fontweight','bold')
                text(0,45,['mean RT: ' num2str(mean(RTs(trialType==2))) ' s'],'color',col1a,'fontsize',16,'fontweight','bold')
                text(0,50,['mean RT: ' num2str(mean(RTs(trialType==3))) ' s'],'color',col1b,'fontsize',16,'fontweight','bold')
                text(0,55,['mean RT: ' num2str(mean(RTs(trialType==4))) ' s'],'color',col2,'fontsize',16,'fontweight','bold')
                hold off
                xlim([-1 2])
                set(gca,'linewidth',2,'fontsize',16)
                
                % for easy trials.
                figure(handleF)
                plotmultipleaxes(2,2,2,0.05,handleF);
                hold on
                % plotting errors
                plot(f,(Sf0dB+Sf0err)','color',col0,'linestyle',':')
                plot(f,(Sf0dB-Sf0err)','color',col0,'linestyle','--')
                plot(f,(Sf1adB+Sf1aerr)','color',col1a,'linestyle',':')
                plot(f,(Sf1adB-Sf1aerr)','color',col1a,'linestyle','--')
                plot(f,(Sf1bdB+Sf1berr)','color',col1b,'linestyle',':')
                plot(f,(Sf1bdB-Sf1berr)','color',col1b,'linestyle','--')
                plot(f,(Sf2dB+Sf2err)','color',col2,'linestyle',':')
                plot(f,(Sf2dB-Sf2err)','color',col2,'linestyle','--')
                % plotting means
                plot(f,Sf0dB,'color',col0,'linewidth',2)
                plot(f,Sf1adB,'color',col1a,'linewidth',2)
                plot(f,Sf1bdB,'color',col1b,'linewidth',2)
                plot(f,Sf2dB,'color',col2,'linewidth',2)
                
                set(gca,'linewidth',2,'fontsize',14,'xscale','log')
                axis tight
                xlim(params.fpass)
                hold off
                
%                 if  sum(isnan(SVec))~=length(SVec)
%                     
%                     figure(handleF)
%                     plotmultipleaxes(4,2,2,0.05,handleF);
%                     hold on
%                     bar(1:size(means,1),means(:,1))
%                     errorbar(1:size(means,1),means(:,1),means(:,2),'linewidth',4,'linestyle','none','color','k')
%                     hold off
%                     set(gca,'fontsize',14,'linewidth',2)
%                     
%                 end
                
                % titling figure with LFP label
                title(deblank(LFPlabels{ch}))
                
                
            case 2
                %% [20180608] 
                spectralType = 'OverallSpectra';
                
                Sf_baseline = nanmean(Sf_baseline,2);
                        
                % collecting mean spectra
                SfdB = nanmean(Sf_encodingWindow./repmat(Sf_baseline,1,nTrials),2);
                                
                % calculating error over trials. 
                Sferr = nanstd(Sf_encodingWindow./repmat(Sf_baseline,1,nTrials),[],2)./sqrt(nTrials);
                
  
                %% saving spectrum data
                fieldStats.Sf_baseline = Sf_baseline; 
                fieldStats.SfdB = SfdB;
                fieldStats.Sferr = Sferr;
                
%                 keyboard
                
                if exist('./Data','dir')
                    try
                        fName = sprintf('./Data/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_SpectraAllTrials',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
                        save([fName '.mat'],'fieldStats')
                    catch
                        mkdir('./Data/LFP/')
                        fName = sprintf('./Data/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_SpectraAllTrials',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
                        save([fName '.mat'],'fieldStats')
                    end
                elseif exist('./Data/LFP/','dir')
                    fName = sprintf('../../Data/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_SpectraAllTrials',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
                    save([fName '.mat'],'fieldStats')
                else
                    fName = sprintf('%s_session_%d_Channel_%s_LFP%s_%sAligned_SpectraAllTrials',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName);
                    save([fName '.mat'],'fieldStats')
                end
                display(sprintf('LFP stats data saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',fName));

                
                %% plotting spectra
                handleF = ch*100;
                brdr = 0.04;
                figure(handleF)

                % for easy trials.
                figure(handleF)
                hold on
                % plotting errors
%                 plot(f,(SfdB+Sferr)','color',col0,'linestyle',':')
%                 plot(f,(SfdB-Sferr)','color',col0,'linestyle','--')
                patch([f fliplr(f)],[SfdB'+Sferr' fliplr(SfdB'-Sferr')],[0.2 0.2 0.2],'edgecolor','none','facealpha',0.3)
                % plotting means
                plot(f,SfdB,'color','k','linewidth',2)

                % deets
                set(gca,'linewidth',2,'fontsize',14,'xscale','log')
                axis tight square
                xlim(params.fpass)
                hold off
                % titling figure with LFP label
                title(deblank(LFPlabels{ch}))
                
        end
        
        fName = sprintf('%s_session_%d_Channel_%s_UMALFP_%sAligned',patientID,sessionNum,deblank(LFPlabels(ch,:)),alignName);
        print(handleF,sprintf('~/data/msit_units/Experiment_II__dlPFC/%s/Figs/spectrograms/%s',patientID,fName),'-dpdf')
        close(handleF)
        
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
%         maximize(handleF)
%         pause(3)
%         if exist('./Figs','dir')
%             try
%                 fName = sprintf('./Figs/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_denoisedWith%s',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignName,denoiseMethod);
%                 saveas(handleF,fName,'pdf')
%                 close(handleF)
%             catch
%                 mkdir('./Figs/LFP/')
%                 fName = sprintf('./Figs/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_denoisedWith%s',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignNamedenoiseMethod);
%                 saveas(handleF,fName,'pdf')
%                 close(handleF)
%             end
%         elseif exist('./Figs/LFP/','dir')
%             fName = sprintf('../../Figs/LFP/%s_session_%d_Channel_%s_LFP%s_%sAligned_denoisedWith%s',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignNamedenoiseMethod);
%             saveas(handleF,fName,'pdf')
%             close(handleF)
%         else
%             fName = sprintf('%s_session_%d_Channel_%s_LFP%s_%sAligned_denoisedWith%s',patientID,sessionNum,deblank(LFPlabels{ch}),spectralType,alignNamedenoiseMethod);
%             saveas(handleF,fName,'pdf')
%             close(handleF)
%         end
%         
%         display(sprintf('figure saved as %s\n\n        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        ',fName));
        
        %
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
        
%         close(gcf)
    end % looping over lfp channels
end % looping over align spots (Stimulus & response)
