function spikeMetrics = unitFeatures(patientID,session,nevFile)
% UNITFEATURES generates summary data for units in a NEV file.
%
%   fName = unitFeatures(patientID, session, nevFile) will generate a
%   pdf that shows unit waveforms and interspike interval histograms for
%   each unit specified in the NEV file.
%
%   The patientID and session arguments are optional.
%
%   The ouput will idicate the location of the pdf file.
%

% author: EHS20160402

close all;

% initializing output variable
spikeMetrics = struct('chan',[],'unit',[],'wfAmplitude',[],'SNR',[]);


%% loading data from NEV file
display('loading unit data...')
% parsing file extension
ext = nevFile(end-3:end);
if strcmp(ext,'.nev')
    NEV = openNEV(nevFile,'read');
elseif strcmp(ext,'.mat')
    load(nevFile);
else % if there isn't an extension specified...
    try
        load([nevFile '.mat']);
    catch
        NEV = openNEV([nevFile '.nev'],'read');
    end
end


%% Looking for task data.
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
timeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% creating neural timing variable
if strcmp(patientID,'CUBF09')
    ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode) double(NEV.Data.Spikes.Unit) double(NEV.Data.Spikes.TimeStamp)];
else
    ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./timeRes)'];
end
ChanUnitTimestamp(ChanUnitTimestamp(:,2)==255,:) = [];

% wtf is the wf
WFs = NEV.Data.Spikes.Waveform./4;
WFtime = linspace(1,(size(WFs,1)-1)./timeRes,size(WFs,1))*1000; % waveform timebase in microseconds.

% chan params
inclChans = unique(ChanUnitTimestamp(:,1));
nChans = length(inclChans);
channelNames = {NEV.ElectrodesInfo.ElectrodeLabel};


%% analyzing over channels and units.
% looping over Channels
for ch = 1:nChans
    % looping over number of units in the AP data
    nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
    for un = 1:nUnits
        
        % getting unit times for the current channel and unit.
        currentUnit = ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un;
        unitTimes = ChanUnitTimestamp(currentUnit,3); % in seconds
        unitWFs = double(WFs(:,currentUnit));
        
        unitWFbar = fliplr(nanmedian(unitWFs'));
        unitWFerr = fliplr(std(double(unitWFs')./sqrt(length(unitWFs)-1)));
        
        
        %% calucalting average firing rate.
        if ~isempty(unitTimes)
            FRbar = length(unitTimes)/(unitTimes(end)-unitTimes(1));
        else
            FRbar = 0;
        end
        
        
        %% plotting waveforms.
        figure(un)
        plotmultipleaxes(1,2,1,0.1,un)
        hold on
        patch([WFtime fliplr(WFtime)],[unitWFbar+unitWFerr fliplr(unitWFbar-unitWFerr)],[0 0 0],'FaceAlpha',0.3)
        plot(WFtime,unitWFbar,'k','linewidth',2)
        
        set(gca,'fontsize',14,'linewidth',2)
        axis square
        xlabel('time (micro s)','fontsize', 16)
        ylabel('mean/std unit waveform (uV)','fontsize', 16)
        
        title(sprintf('channel %s unit %d average firing rate = %d spikes/s',channelNames{ch}(1:5),un,FRbar))
        
        
        %% plotting ISI histogram.
        figure(un)
        % meats
        plotmultipleaxes(2,2,1,0.1,un)
        histogram(diff(unitTimes),0:0.005:0.5,'DisplayStyle','stairs')

        % deets
        set(gca,'fontsize',14,'linewidth',2)
        axis square
        xlabel('time (s)','fontsize', 16)
        ylabel('AP count','fontsize', 16)
        title(sprintf('total spikes = %d',length(unitTimes)))
        
        
        %% saving stats and figs.
        display('saving stats and figs to Elliot"s dropbox. Thanks!')
        dbPath = '/media/user1/data4TB/msit_units/unitData';
        % output data
        spikeMetrics.chan = spikeMetrics.chan;
        spikeMetrics.unit = spikeMetrics.unit;
        spikeMetrics.wfAmplitude = max(unitWFbar)-min(unitWFbar);
        spikeMetrics.SNR = NEV.ElectrodesInfo(inclChans(ch)).LowThreshold./4;
        fName = sprintf('%s/%s_session_%d_Channel_%d_Unit_%d_%s_waveformStats.mat',dbPath,patientID,session,inclChans(ch),un,deblank(channelNames{ch}'));
        save(fName,'spikeMetrics')
        
%         % saving figures.
%         fName = sprintf('%s/%s_session_%d_Channel_%d_Unit_%d_%s_waveforms',dbPath,patientID,session,inclChans(ch),un,deblank(channelNames{ch}'));
%         saveas(un,fName, 'pdf')
        close(un)
            
            
    end % looping over units for each channel and align spot
end % looping over channels for each align spot.

