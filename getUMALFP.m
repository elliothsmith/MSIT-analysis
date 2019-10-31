% grabs LFP recorded on CUCX2's microelectrode array
clear
clc

%% load data
for fl = 1:3
    switch fl
        case 1
            dataFile = './Experiment_II__dlPFC/CUCX2/Data/20160802-111303-001.ns5';
        case 2
            dataFile = './Experiment_II__dlPFC/CUCX2/Data/20160803-115008-001.ns5';
        case 3
            dataFile = './Experiment_II__dlPFC/CUCX2/Data/20160807-140905-008.ns5';            
    end
    fprintf('\nfile: %d',fl)
    clearvars -except fl dataFile 
    
    %% loading data. 
    openNSx(dataFile);
    
    
    %% for filtering and resampling each channel
    % first resampling each channel. 
    nChans = 96;
    b = fir1(90,100/(3e4/2));
    for ch = 1:nChans
        fprintf('\n    channel: %d',ch)
        dfData(ch,:) = filtfilt(b,1,double(NS5.Data{2}(ch,:)));
    end
    clear NS5
    
    % resampling down-filtered signal. 
    dsData = nanmean(resample(dfData',500,3e4)');
    
    
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
    plot(Fuma,smooth(10*log10(Suma),50),'color','r')
    plot(Fecog,smooth(10*log10(Secog),50),'color','b')
    hold off
    axis square
    
    % time-domain correlation
    subplot(2,1,2)
    [rho,p] = corrcoef(dsData',dsECoG(1:length(dsData))');
    hold on
    plot(dsData,'color','r')
    plot(dsECoG,'color','b')
    text(1e3,1e3,sprintf('corrcoef: %d, pval: %d',rho(1,2),p(1,2))) 
    hold off
    
    maximize(gcf)
    saveas(gcf,sprintf('correlationBetweenUMAandECoGLFPs_%d.fig',fl))
    
end