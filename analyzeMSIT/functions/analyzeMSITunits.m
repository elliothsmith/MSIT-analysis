function [neuronStats,Rates] = analyzeMSITunits(patientID,sessionNum,nevFile)
%ANALYZEMSITUNITS plots PSTHs for MSIT units
%
%   [neuronStats] = analyzeMSITunits(patientID,sessionNumber,nevFile) will
%       plot PSTHs for sorted units recorded in a blackrock NEV file or its
%       associated .mat and return statistics.
%
%   neuronStats contains
%

% Author: EHS 20160712
% VersionControl: https://github.com/elliothsmith/MSIT-analysis


%% [20160831] initializing ouptut variable
Rates = struct();


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
if strcmp(patientID,'CUBF09')
    ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode) double(NEV.Data.Spikes.Unit) double(NEV.Data.Spikes.TimeStamp)];
else
    ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];
end
ChanUnitTimestamp(ChanUnitTimestamp(:,2)==255,:) = [];


inclChans = unique(ChanUnitTimestamp(:,1));
nChans = length(inclChans);


%% aligning AP data.
display('Aligning AP data on stimulus and response.');
% for loop to save multiple epochs
for aS = 1:2
    % which alignment spot
    switch aS
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            %% timing (seconds)
            pre = 2;
            post = 3;
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>=100 & trigs<=105);
            %% timing (seconds)
            pre = 3;
            post = 3;
    end
    
    
    % looping over Channels
    for ch = 1:nChans
        % looping over number of units in the AP data
        nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
        for un = 1:nUnits
            % getting unit times for the current channel and unit.
            unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds
            % loooping over trials
            for tt = 1:nTrials
                
                %% putting the data in a structure
                data(aS).channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post) - repmat(trialStarts(tt)-pre,length(unitTimes(unitTimes>trialStarts(tt)-pre & unitTimes<trialStarts(tt)+post)),1);
                
            end
        end
    end
    
    % save data structure?
    saveFlag = 0;
    if saveFlag
        display(['saving spike time structure for ' alignName '...'])
        save([patientID '_session' num2str(sessionNum) '_spikeTimeStruct_alignedon' alignName '.mat'],'data','pre','post','alignName')
    end
    
end


%% parsing behavior
condition = trigs(trigs>=1 & trigs<=27);
easyTrials = condition>=1 & condition<=3;
simonTrials = condition>=4 & condition<=21;
eriksenTrials = (condition>=4 & condition<=15) | (condition>=22 & condition<=27);


%% Rasters and PSTHs
for aS2 = 1:length(data)
    % which alignment spot
    switch aS2
        case 1
            alignName = 'Cue';
            trialStarts =  trigTimes(trigs>=1 & trigs<28);
            %% timing (seconds)
            pre = 2;
            post = 3;
        case 2
            alignName = 'Response';
            trialStarts =  trigTimes(trigs>=100 & trigs<=105);
            %% timing (seconds)
            pre = 3;
            post = 3;
    end
    
    % looping over units
    for ch = 1:size(data(aS2).channel,2)
        for un = 1:size(data(aS2).channel(ch).unit,2)
            
            display('plotting raster over conflict... grab a cocktail, this may take a while.')
            %% plotting rasters and PSTHs
            figure(aS2)
            ah_ras = plotmultipleaxes(1,1,2,0.08,aS2);
            hold on
            for tt = 1:nTrials
                % changing raster color based on trial type
                if easyTrials(tt)
                    rasCol = col0;
                elseif simonTrials(tt)
                    rasCol = col1a;
                elseif eriksenTrials(tt)
                    rasCol = col1b;
                end
                
                % plotting rasters for conflict (in the least efficient way possible)
                for sp = 1:size(data(aS2).channel(ch).unit(un).trial(tt).times,1)
                    try
                        line([data(aS2).channel(ch).unit(un).trial(tt).times(sp)-pre data(aS2).channel(ch).unit(un).trial(tt).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', rasCol)
                    catch
                        line([data(aS2).channel(ch).unit(un).trial(tt).times(sp)-pre data(aS2).channel(ch).unit(un).trial(tt).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', 'k')
                    end
                end
                
                
                % stimulus timing lines
                line([0 0], [0 nTrials],'linestyle', '--', 'color', 'k')
                % raster plot details
                xlim([-(pre-1) (post-1)])
                ylim([0 nTrials])
                str = sprintf('patient %s, Channel %d, Unit %d; aligned on %s',patientID ,ch ,un ,alignName);
                title(str,'fontsize',18);
                ylabel('Trials','fontsize', 16)
                set(gca, 'linewidth', 2, 'fontsize', 16);
                
                
            end
            hold off
            
            
            %% PSTHs
            % calculating psths
            kernelWidth = 50  ./1000;
            [Reasy,t,Eeasy] = psth(data(aS2).channel(ch).unit(un).trial(easyTrials), kernelWidth, 'n', [0 pre+post]);            
            
            try
                [Rsimon,t,Esimon] = psth(data(aS2).channel(ch).unit(un).trial(simonTrials), kernelWidth, 'n', [0 pre+post]);
            catch
                display('no spatial conflict')
            end
            
            try
                [Reriksen,t,Eeriksen] = psth(data(aS2).channel(ch).unit(un).trial(eriksenTrials), kernelWidth, 'n', [0 pre+post]);
            catch
                display('no distracter conflict')
            end
            
            %         % overall PSTH
            %         [R,t,E] = psth(data.channel(ch).unit(un).trial, kernelWidth, 'n', [-pre post]);
            tsec = t-repmat(pre,1,length(t));
            
            
            %% generate statistics for each Neuron.
            neuronStats = [];
%             [pEasy,Heasy] = ranksum(Reasy(tsec>=-1 & tsec<=-0.5),Reasy(tsec>=0 & tsec<=1.5));
%             if isequal(aS,1)
%                 neuronStats.Cue{ch,un} = table(pEasy,pHard,'VariableNames',{'pEasyvsBaseline','pHardvsBaseline'});
%             elseif isequal(aS,2)
%                 neuronStats.Response{ch,un} = table(pEasy,pHard,'VariableNames',{'pEasy','pHard'});
%             end
            
            %% [2016083] do GLM here
            
            
            %% PSTH plots
            ah_psth = plotmultipleaxes(2,1,2,0.08,aS2);
            hold on
            try
                % plotting med psth
                patch([tsec fliplr(tsec)],[Rsimon+Esimon fliplr(Rsimon-Esimon)], col1a,'edgecolor','none','facealpha',0.5)
                plot(tsec,Rsimon,'color',col1a,'linewidth',2)
            catch
                display('no spatial conflict')
            end
            
            try
                % plotting med psth
                patch([tsec fliplr(tsec)],[Reriksen+Eeriksen fliplr(Reriksen-Eeriksen)], col1b,'edgecolor','none','facealpha',0.5)
                plot(tsec,Reriksen,'color',col1b,'linewidth',2)
            catch
                display('no spatial conflict')
            end
            
            % plotting easy psth
            patch([tsec fliplr(tsec)],[Reasy+Eeasy fliplr(Reasy-Eeasy)], col0,'edgecolor','none','facealpha',0.5)
            plot(tsec,Reasy,'color',col0,'linewidth',2)
                        
            % stimulus timing lines
            line([0 0], [0 nTrials],'linestyle', '--', 'color', 'k')
            % PSTH plot details
            xlim([-(pre-1) (post-1)])
            ylim([0 max(Reriksen+Eeriksen)+3])
            hold off
            xlabel('Time (seconds)', 'fontsize', 16);
            ylabel('Firing Rate (spikes/second)', 'fontsize', 16);
            set(gca, 'linewidth', 2, 'fontsize', 16);
            
            
            %% [20160831] organizing output variable
            Rates.tsec = tsec;
            Rates.conflict.none.rate = Reasy;
            Rates.conflict.none.error = Eeasy;
            Rates.conflict.simon.rate = Rsimon;
            Rates.conflict.simon.error = Esimon;
            Rates.conflict.flanker.rate = Reriksen;
            Rates.conflict.flanker.error = Eeriksen;
            
            keyboard
            
            %% saving figures.
            figFlag = 1;
            if figFlag
                if exist('./Figs','dir')
                    if exist('./Figs/FiringRate/','dir')
                        fName = sprintf('./Figs/FiringRate/%s_session_%d_Channel_%d_Unit_%d_SimonEriksenConflict_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
                        saveas(aS2,fName, 'pdf')
                        close(aS2)
                    else
                        mkdir('./Figs/FiringRate/')
                        fName = sprintf('./Figs/FiringRate/%s_session_%d_Channel_%d_Unit_%d_SimonEriksenConflict_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
                        saveas(aS2,fName, 'pdf')
                        close(aS2)
                    end
                else
                    try
                        fName = sprintf('./Figs/%s_session_%d_Channel_%d_Unit_%d_SimonEriksenConflict_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
                        saveas(aS2,fName, 'pdf')
                        close(aS2)
                    catch
                        fName = sprintf('%s_session_%d_Channel_%d_Unit_%d_SimonEriksenConflict_%saligned',patientID,sessionNum,inclChans(ch),un,alignName);
                        saveas(aS2,fName, 'pdf')
                        close(aS2)
                    end
                end
            end
            
            
%             %% saving stats.
%             if exist('./Data','dir')
%                 fName = sprintf('./Data/%s_session_%d_ConflictStats_alignedOn_%s',patientID,sessionNum,alignName);
%                 save([fName '.mat'],'neuronStats')
%             elseif ~exist('./Data','dir')
%                 mkdir('./Data/')
%                 fName = sprintf('./Data/%s_session_%d_ConflictStats_alignedOn_%s',patientID,sessionNum,alignName);
%                 save([fName '.mat'],'neuronStats')
%             else
%                 fName = sprintf('%s_session_%d_ConflictStats_alignedOn_%s',patientID,sessionNum,alignName);
%                 save([fName '.mat'],'neuronStats')
%             end
            
            
        end % looping over units for each channel and align spot
    end % looping over channels for each align spot.
end % looping over align spots (Stimulus & response)

