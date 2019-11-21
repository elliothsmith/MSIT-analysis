% plot unit scatter
close all
 
%% loading GLM unit classification data.
msitUnitsPath = '/media/user1/data4TB/';
load([msitUnitsPath 'msit_units/acc_dlpfc_units_results.mat'])
[acc_dlpfc_units,acc_units,dlpfc_units] = parseMSITGLMresults();


%% finding the right units and plotting.
unitDir = '/media/user1/data4TB/msit_units/unitDataACC';
dirList = dir(unitDir);
dirList = dirList(3:end);

unitIdcs = find(acc_units);


%% looping over data points.
for p = 1:136
    % loading data
    load(fullfile(unitDir,dirList(unitIdcs(p)).name))
    
    % plotting data for dACC
    figure(1)
    subplot(1,2,1)
    hold on
    scatter(abs(spikeMetrics.SNR),spikeMetrics.wfAmplitude,3,[0 0 0],'filled');
    axis square
    xlim([10 100])
    ylim([10 100])
    hold off
    
end
xlabel('channel threshold (uV)')
ylabel('unit amplitude (uV)')
title('dACC')

%% finding the right units and plotting.
unitDir = '/media/user1/data4TB/msit_units/unitDataPFC';
dirList = dir(unitDir);
dirList = dirList(3:end);


%% looping over data points.
for p = 1:136
    % loading data
    load(fullfile(unitDir,dirList(unitIdcs(p)).name))
    
    % plotting data for dACC
    figure(1)
    subplot(1,2,2)
    hold on
    scatter(abs(spikeMetrics.SNR),spikeMetrics.wfAmplitude,3,[0 0 0],'filled');
    axis square
    xlim([10 100])
    ylim([10 100])
    hold off

end
xlabel('channel threshold (uV)')
ylabel('unit amplitude (uV)')
title('dlPFC')

saveas(1,'~/Dropbox/ACClPFCunitMetrics.pdf')


