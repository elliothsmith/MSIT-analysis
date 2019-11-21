function sessionStats = analyzeMSIToverSessions(patientID,startSession,varargin)
% ANALYZEMSITOVERSESSIONS MSIT analysis over sessions for a single subject.
%   Finds and enters directory called [PatientID].
%   Looks for data in subdirectory ./Data and analyzes the data in that
%   directory across behavioral sessions.
%
%   The code will start analyzing from [startSession] after the first
%   session. i.e. startSession=1 => starts from 2nd session.
%
%   The only required input argument is the patient ID. (e.g. 'CUBF09')
%
%   Other input arguments are a series of strings corresponding to
%   analyses one would like to perform and the functions in the toolbox.
%
%   The potential inputs are:
%       1) behavior
%       2) units
%       3) fields
%       4) coherence
%       5) fieldCMI
%       6) waves
%       7) aptep
%       8) buttonFBtCode ('button and feedback coherence control')
%
%       for example: analyzeMSIToverSessions('CUBF10',0,'behavior','units')
%
%   This function thus acts as a wrapper for the other functions in the
%   MSIT-analysis toolbox, and again works great if data is organized into
%   patient directories with sub-directories called '/patientID/Figs/' (for
%   figures) and '/patientID/Data/' (for data).
%

% Author: Elliot H Smith
% Version Control: https://github.com/elliothsmith/MSIT-analysis


%% [20160711] getting input arguments.
argList = varargin';

sessionStats = struct();

%% TODO: get patient Directory.
% if isequal(exist(['./' patientID],'dir'),7)
%     display('patient directory exists in cwd.')
%     patientDir = './';
% else
%     patientDir = uigetdir('~/Data_B/msit_units/','Select directory with patient sub-directories');
%     cd(patientDir)
% end

if strcmp(patientID,'CUCX2')
    patientDir = '/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/'
else
    patientDir = '/media/user1/data4TB/data/msit_units/Experiment_I__ACC/'
end
cd(patientDir)

%% running the functions from the input arguments.
if ~isequal(exist(patientID,'dir'),7)
    display('no patient directory found')
else
    % entering patient directory
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
        dirList = dir(fullfile(fullPath,'Data/','*NEV.mat'));
        display(['found ' num2str(length(dirList)) ' .mattified NEV files!'])
    end
    
    % looping over sessions
    for  ss = 1+startSession:length(dirList)
        nevFile = fullfile(fullPath,'Data',dirList(ss).name);
        
        % analyzing behavior and saving figures and stats
        if cell2mat(strfind(argList,'behavior'))
            [conflictStats{ss},FBstats{ss},RTs{ss},trialType{ss},errors{ss}]...
                = analyzeMSITbehavior(patientID,ss,nevFile);
            
            % saving behavioral statistics.
            save(['./Data/' patientID '_behavioralStatistics_OverSessions.mat'],'conflictStats','FBstats','RTs','trialType','errors')
        end
        
        % analyzing neurons and saving figures and stats
        if cell2mat(strfind(argList,'units'))
            [neuralStats{ss},Rates{ss}] = analyzeMSITunits (patientID,ss,nevFile);
            
            % saving firing rate statistics.
            save(['./Data/' patientID '_firingRateStatistics_OverSessions.mat'],'neuralStats','Rates')
        end
        
        % analyzing LFP and saving figures and stats
        if cell2mat(strfind(argList,'fields'))
            % for the NatNeuro appeal
            appeal = true;
            if appeal
                if strcmp(patientID,'CUCX2')
                    [fieldStats] = analyzeMSITfieldsU_appeal(patientID,ss,nevFile,1);
                else
                    [fieldStats] = analyzeMSITfields_appeal(patientID,ss,nevFile,1);
%                     [LFPlabels] = analyzeMSITfields_LABELS(patientID,ss,nevFile,1);
                end
            else
                % original analyses.
                if strcmp(patientID,'CUCX2')
                    [fieldStats] = analyzeMSITfieldsU(patientID,ss,nevFile,1);
                else
                    [fieldStats] = analyzeMSITfields (patientID,ss,nevFile,1);
                end
            end
            % saving LFP statistics.
            save(['./Data/' patientID '_LFPStatistics_OverSessions.mat'],'fieldStats')
%             save(['./Data/' patientID '_LFPtrodeLabels.mat'],'LFPlabels')
        end
        
        % analyzing LFP WAVES and saving figures and stats
        if cell2mat(strfind(argList,'waves'))
            if strcmp(patientID,'CUCX2')
                analyzeMSITwavesU (patientID,ss,nevFile);
            else
                analyzeMSITwaves (patientID,ss,nevFile);
            end
            
            %             % saving LFP statistics.
            %             save(['./Data/' patientID '_LFPStatistics_OverSessions.mat'],'fieldStats')
        end
        
        if cell2mat(strfind(argList,'populationCoherence'))
            analyzeMSITpopulationCoherenceU (patientID,ss,nevFile);
        end
        
        
        % analyzing AP-LFP coherence and saving figures and stats
        if cell2mat(strfind(argList,'coherence'))
            if strcmp(patientID,'CUCX2')
                analyzeMSITcoherenceU (patientID,ss,nevFile);
            else
                analyzeMSITcoherence (patientID,ss,nevFile);
            end
            
            % saving LFP statistics.
            %             save(['./Data/' patientID '_coherenceStatistics_OverSessions.mat'],'coherenceStats')
        end
        
        % analyzing LFP CMI and saving figures and stats
        if cell2mat(strfind(argList,'fieldCMI'))
            analyzeMSITfieldCMI (patientID,ss,nevFile);
            
            %             % saving LFP statistics.
            %             save(['./Data/' patientID '_CMIStatistics_OverSessions.mat'],'CMIStats')
        end
        
        % analyzing stLFP and saving figures and stats
        if cell2mat(strfind(argList,'aptep'))
            if strcmp(patientID,'CUCX2')
                analyzeMSITaptepU (patientID,ss,nevFile);
            else
                analyzeMSITaptep (patientID,ss,nevFile);
            end
            
            % saving LFP statistics.
            %             save(['./Data/' patientID '_coherenceStatistics_OverSessions.mat'],'coherenceStats')
        end
        
        % analyzing AP-LFP coherence and saving figures and stats
        if cell2mat(strfind(argList,'buttonFBtCode'))
            if strcmp(patientID,'CUCX2')
                display('...')
                %                 analyzeMSITcoherenceU (patientID,ss,nevFile);
            else
                analyzeMSITcoherence_buttonAndFB(patientID,ss,nevFile);
            end
            
            % saving LFP statistics.
            %             save(['./Data/' patientID '_coherenceStatistics_OverSessions.mat'],'coherenceStats')
        end
        
        % analyzing LFP-LFP coherence and saving figures and stats
        if cell2mat(strfind(argList,'fieldFieldCoherence'))
            if strcmp(patientID,'CUCX2')
                display('...')
                %                 analyzeMSITcoherenceU (patientID,ss,nevFile);
            else
                analyzeMSITfieldFieldCoherence(patientID,ss,nevFile);
            end
            
            % saving LFP statistics.
            %             save(['./Data/' patientID '_coherenceStatistics_OverSessions.mat'],'coherenceStats')
        end
        
        
        % analyzing neurons and saving figures and stats
        if cell2mat(strfind(argList,'unitFeatures'))
            if ss==1
                wfAmps = [];
            end
            spikeMetrics = unitFeatures(patientID,ss,nevFile);
            
            wfAmps = cat(1,wfAmps,spikeMetrics.wfAmplitude);
        end
        
        %% TODO: saving statistics over sessions. === Do this in each section. makes more sense.
        %         saveFlag = 1;
        %         if saveFlag
        %
        %         end
        
    end
    % goping back
    cd(patientDir)
end

%% TODO:: create summary report.


if cell2mat(strfind(argList,'unitFeatures'))
    sessionStats(ss).waveformAmplitudes = wfAmps;
end



end

