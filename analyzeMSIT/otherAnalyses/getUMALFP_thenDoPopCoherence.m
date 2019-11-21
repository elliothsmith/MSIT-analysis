% grabs LFP recorded on CUCX2's microelectrode array
clear
clc

% two steps. 
getLFP = false;
calcPopCoh = true;


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
        
        
        %% for down-sampling the mean across electrodes.
        %     % taking the mean across electrodes
        %     meanData = median(double(NS5.Data{2}(1:92,:)));
        %     clear NS5
        %
        %     % downsampling data to 500 Hz.
        %     dsData = resample(meanData,500,3e4);
        
        
        %     %% saving data
        save([dataFile(1:end-4) '_downsampledMeanUMALFP.mat'],'dsData','-v7.3')
        
        
        % loading ECoG data.
        openNSx([dataFile(1:end-1) '3'])
        
        % downsampling EcoG channel
        dsECoG = resample(double(NS3.Data{2}(56,:)),500,2e3);
        
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
end

% calculating coherence
if calcPopCoh
    %% load data
    for fl = 3
        switch fl
            case 1
                nevFile = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/*.mat';
            case 2
                nevFile = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/*.mat';
            case 3
                nevFile = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/*.mat';
        end
        disp(fl)
        
        % do coherence calculations.
        analyzeMSITpopulationCoherenceU_UMALFP('CUCX2',fl,nevFile)
        
        
    end
end



