% grabs LFP recorded on CUCX2's microelectrode array
clear
clc
close all

% 1) derive the analogous LFP for the UMA
getLFP = false;

% 2) comparisons among recording modalitites.
compareUMAandECoG = false;
compareECoGandsEEG = false;
withinPtDiffs = false; % previous two need to be true for this one to work...
broadbandLFPcorrelations = true;

% 3) calculate population coherence with UMA LFP.
calcPopCoh = false;


%% doing LFP preprocessing.
if getLFP
    %% load data
    for fl = 1:3
        switch fl
            case 1
                dataFile = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/******02-111303-001.ns5';
            case 2
                dataFile = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/******03-115008-001.ns5';
            case 3
                dataFile = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/******07-140905-008.ns5';
        end
        fprintf('\nfile: %d',fl)
        clearvars -except fl dataFile getLFP calcPopCoh
        
        
        %% loading data.
        openNSx(dataFile);
        
        
        %% for filtering and resampling each channel
        % first resampling each channel.
        nChans = 96;
        b = fir1(90,50/(3e4/2));
        for ch = 1:nChans
            fprintf('\n    channel: %d',ch)
            dfData(ch,:) = filtfilt(b,1,double(NS5.Data{2}(ch,:)));
        end
        clear NS5
        
        %         % visualize data. There may be some bad channels in there.
        %         figure
        %         imagesc(dfData,[-1e3 1e3])
        %         keyboard
        
        
        % resampling, removing common mode, and taking the mean of the down-filtered signal.
        dsData = nanmean((resample(dfData(11:80,:)',500,3e4)'));
        
        
    end
end



if compareUMAandECoG
    %% load data
    for fl = 3
        switch fl
            case 1
                dataFile = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/******02-111303-001.ns5';
            case 2
                dataFile = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/******03-115008-001.ns5';
            case 3
                lastTrial = 289;
                dataFile = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/******07-140905-008.ns5';
        end
        
        % loading the preprocessed UMA LFP
        load([dataFile(1:end-4) '_downsampledMeanUMALFP.mat'])
        
        % loading and downsampling ECoG data.
        openNSx([dataFile(1:end-1) '3'])
        tmp1 = remove1stPC(double(NS3.Data{2}));
        dsECoG = resample(tmp1(56,:),500,2e3);
        
        % which comparisons to examine
        overallComparison = false;
        trialComparison = true;
        
        if overallComparison
            % plotting LFP comparison
            figure
            % overall frequecny domain comparisons
            % [20190122] this is comparing the overall frequencies, but we could
            % (paired) test the spectra across trials for statistics.
            subplot(2,2,1)
            [Suma,Fuma] = periodogram(dsData,hann(length(dsData)),1e4,500);
            [Secog,Fecog] = periodogram(dsECoG,hann(length(dsECoG)),1e4,500);
            hold on
            plot(Fuma,zscore(smooth(10*log10(Suma),50)),'color','r')
            plot(Fecog,zscore(smooth(10*log10(Secog),50)),'color','b')
            hold off
            axis square
            
            % time-domain correlation
            subplot(2,1,2)
            [rho,p] = corrcoef(dsData',dsECoG(1:length(dsData))');
            hold on
            plot(zscore(dsECoG),'color','b')
            plot(zscore(dsData),'color','r')
            text(1e3,1e3,sprintf('corrcoef: %d, pval: %d',rho(1,2),p(1,2)))
            hold off
            
            maximize(gcf)
            saveas(gcf,sprintf('correlationBetweenUMAandECoGLFPs_%d.fig',fl))
        end
        if trialComparison
            % get nev file
            nevList = dir([dataFile(1:end-4) '*NEV*'])
            load(fullfile(nevList.folder,nevList.name))
            
            % get markers
            trigs = NEV.Data.SerialDigitalIO.UnparsedData;
            trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
            cueTimes =  trigTimes(trigs>=1 & trigs<28);
            
            % timing params
            pre = 3;
            post = 5;
            Fs = 500;
            params.fpass = [1 50];
            
            % epoch data
            for tt = 1:lastTrial
                UMALFP(:,tt) = dsData(floor((cueTimes(tt)-pre)*Fs):floor((cueTimes(tt)+post)*Fs));
                ECoGLFP(:,tt) = dsECoG(floor((cueTimes(tt)-pre)*Fs):floor((cueTimes(tt)+post)*Fs));
                
                [wUMA,periodUMA,~] = basewaveERP(UMALFP(:,tt),Fs,params.fpass(1),params.fpass(2),6,0);
                [wECoG,periodECoG,~] = basewaveERP(ECoGLFP(:,tt),Fs,params.fpass(1),params.fpass(2),6,0);
                
                % power
                Sft_UMA(:,:,tt) = abs(wUMA);
                Sft_ECoG(:,:,tt) = abs(wECoG);
            end
            
            % determining spectral scales in Hz
            freqsUMA = periodUMA.^-1;
            freqsECoG = periodECoG.^-1;
            
            % timing stuff
            tSec = linspace(-pre,post,500*(pre+post)+1);
            bP = [-.9 -.4];
            
            baselineNormalize = true
            if baselineNormalize
                % baseline normalize UMA
                tmp = Sft_UMA./repmat(mean(mean(Sft_UMA(:,tSec>bP(1) & tSec<bP(2),:),2),3),1,size(Sft_UMA,2),size(Sft_UMA,3));
                Sft_UMA = tmp;
                clear tmp
                % baseline normalize ECoG
                tmp = Sft_ECoG./repmat(mean(mean(Sft_ECoG(:,tSec>bP(1) & tSec<bP(2),:),2),3),1,size(Sft_ECoG,2),size(Sft_ECoG,3));
                Sft_ECoG = tmp;
                clear tmp
            end
            
            % stats!!!!!
            timeIdcs = (tSec>0 & tSec<2);
            foundIdcs = find(timeIdcs);
            for ts = 1:sum(timeIdcs)
                [hUMAECoG,pUMAECoG] = ttest(UMALFP(foundIdcs(ts),:),ECoGLFP(foundIdcs(ts),:));
            end
            adjustedP = 1-((1-0.05).^(1/sum(timeIdcs)))
            
        end
    end
end



if compareECoGandsEEG
    %% load data
    lastTrial = 289;
    dataFile = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/******07-140905-008.ns5';
    
    % loading and resampling ECoG data.
    openNSx([dataFile(1:end-1) '3'])
    tmp1 = remove1stPC(double(NS3.Data{2}));
    dsECoG = resample(tmp1(56,:),500,2e3);
    
    % ECoG markers
    nevList = dir([dataFile(1:end-4) '*NEV*'])
    load(fullfile(nevList.folder,nevList.name))
    
    % get markers
    ecog_trigs = NEV.Data.SerialDigitalIO.UnparsedData;
    ecog_trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
    ecog_cueTimes =  ecog_trigTimes(ecog_trigs>=1 & ecog_trigs<28);
    
    % loading and resampling sEEG data.
    pt = 2;
    switch pt
        case 1
            pt = 'CUBF09';
            elec = 15;
        case 2
            pt = 'CUBF17';
            elec = 24;
    end
    tmp1 = remove1stPC(double(NS3.Data));
    dssEEG = resample(tmp1(elec,:),500,2e3);
    
    % sEEG markers
    seeg_trigs = NEV.Data.SerialDigitalIO.UnparsedData;
    seeg_trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
    seeg_trigs(95) = 104;
    seeg_cueTimes =  seeg_trigTimes(seeg_trigs>=1 & seeg_trigs<28);
    %[~,trigTimes] = removePracticeTriggers(sEEG_trigs,trigTimes);
    
    % which comparisons to examine
    trialComparison = true;
    if trialComparison
        % timing params
        pre = 3;
        post = 5;
        Fs = 500;
        params.fpass = [1 50];
        
        % epoch data
        for tt = 1:lastTrial
            seegLFP(:,tt) = dssEEG(floor((seeg_cueTimes(tt)-pre)*Fs):floor((seeg_cueTimes(tt)+post)*Fs));
            ECoGLFP(:,tt) = dsECoG(floor((ecog_cueTimes(tt)-pre)*Fs):floor((ecog_cueTimes(tt)+post)*Fs));
            
            [wseeg,periodsEEG,~] = basewaveERP(seegLFP(:,tt),Fs,params.fpass(1),params.fpass(2),6,0);
            [wECoG,periodECoG,~] = basewaveERP(ECoGLFP(:,tt),Fs,params.fpass(1),params.fpass(2),6,0);
            
            % power
            Sft_sEEG(:,:,tt) = abs(wseeg);
            Sft_ECoG(:,:,tt) = abs(wECoG);
        end
        
        % determining spectral scales in Hz
        freqssEEG = periodsEEG.^-1;
        freqsECoG = periodECoG.^-1;
        
        % timing stuff
        tSec = linspace(-pre,post,500*(pre+post)+1);
        bP = [-.9 -.4];
        
        baselineNormalize = true
        if baselineNormalize
            % baseline normalize UMA
            tmp = Sft_sEEG./repmat(mean(mean(Sft_sEEG(:,tSec>bP(1) & tSec<bP(2),:),2),3),1,size(Sft_sEEG,2),size(Sft_sEEG,3));
            Sft_sEEG = tmp;
            clear tmp
            % baseline normalize ECoG
            tmp = Sft_ECoG./repmat(mean(mean(Sft_ECoG(:,tSec>bP(1) & tSec<bP(2),:),2),3),1,size(Sft_ECoG,2),size(Sft_ECoG,3));
            Sft_ECoG = tmp;
            clear tmp
        end
        
        % stats!!!!!
        timeIdcs = (tSec>0 & tSec<2);
        foundIdcs = find(timeIdcs);
        for ts = 1:sum(timeIdcs)
            [hsEEGECoG,psEEGECoG] = ttest(seegLFP(foundIdcs(ts),:),ECoGLFP(foundIdcs(ts),:));
        end
        adjustedP = 1-((1-0.05).^(1/sum(timeIdcs)))
        
        
    end
end

if withinPtDiffs
    cLim = [0.8 1.2];
    xLims = [-1 3];
    
    figure
    subplot(3,3,1)
    imagesc(tSec,freqsUMA,mean(Sft_UMA,3),cLim+[0 0.1])
    colorbar
    axis xy square
    title('UMA LFP')
    ylabel('LFP frequency (Hz)')
    xlabel('time relative to cue (s)')
    xlim(xLims)
    
    subplot(3,3,2)
    imagesc(tSec,freqsECoG,mean(Sft_ECoG,3),cLim)
    colorbar
    axis xy square
    title('ECoG LFP')
    ylabel('LFP frequency (Hz)')
    xlabel('time relative to cue (s)')
    xlim(xLims)
    
    subplot(3,3,3)
    imagesc(tSec,freqssEEG,mean(Sft_sEEG,3),cLim)
    colorbar
    axis xy square
    title('sEEG LFP')
    ylabel('LFP frequency (Hz)')
    xlabel('time relative to cue (s)')
    xlim(xLims)
    
    subplot(3,3,5)
    imagesc(tSec,freqsECoG,mean(Sft_UMA,3)-mean(Sft_ECoG,3),[-0.2 0.2])
    colorbar
    axis xy square
    title('UMA - ECoG LFP')
    ylabel('LFP frequency (Hz)')
    xlabel('time relative to cue (s)')
    xlim(xLims)
    
    subplot(3,3,6)
    imagesc(tSec,freqssEEG,mean(Sft_UMA,3)-mean(Sft_sEEG,3),[-0.2 0.2])
    colorbar
    axis xy square
    title('UMA - sEEG LFP')
    ylabel('LFP frequency (Hz)')
    xlabel('time relative to cue (s)')
    xlim(xLims)
    
    
    subplot(3,3,9)
    imagesc(tSec,freqssEEG,mean(Sft_ECoG,3)-mean(Sft_sEEG,3),[-0.2 0.2])
    colorbar
    axis xy square
    title('ECoG - sEEG LFP')
    ylabel('LFP frequency (Hz)')
    xlabel('time relative to cue (s)')
    xlim(xLims)
    
    print(gcf,'~/Dropbox/comparingAllLFPdlPFCspectrograms','-dpdf','-fillpage')
end


