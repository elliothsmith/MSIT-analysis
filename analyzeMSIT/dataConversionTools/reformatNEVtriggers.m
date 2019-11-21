function [reconstructedNEV] = reformatNEVtriggers(nevFile,resultMatFile)
% REFORMATNEVTRIGGERS reformats any corrupted triggers in a nev file by 
%   replacing them with the codes that were saved in the .mat file from
%   psych toolbox.


% Author: EHS 20160812
% VersionControl: https://github.com/elliothsmith/MSIT-analysis


%% loading data from NEV file
display('loading data...')
% [pathstr, name, ext] = fileparts(nevFile);

ext = nevFile(end-3:end);
if strcmp(ext,'.nev')
    NEV = openNEV(nevFile,'read');
elseif strcmp(ext,'.mat')
    load(nevFile);
end

% loading behavioral data from psychtoolbox output.
load(resultMatFile)


%% organizing behavioral markers
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;


%% [20160824] removing practice trial triggers and updating number of trials.
[trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
nTrials = sum(trigs==90);


%% most frequently occuring value
% finding starting value
newStartValue = mode(trigs);
% updating starting values
trigs(trigs==newStartValue) = 90;
trigs(trigs==newStartValue+1) = 91;
trigs(trigs==newStartValue+2) = 92;
% updating condition values
startIdcs = find(trigs==92);
conditionIdcs = startIdcs+repmat(1,length(startIdcs),1);
trigs(conditionIdcs) = resultMat(1:length(conditionIdcs),6);
% updating reponse values
responseIdcs = startIdcs+repmat(2,length(startIdcs),1);
trigs(responseIdcs) = resultMat(1:length(conditionIdcs),7) + repmat(100,length(conditionIdcs),1);


if strcmp(nevFile,'/home/elliot/Data_B/msit_units/CUCX2/Data/preformattedNEVs/****0803-115008-001-01_sortedEHS****0811.nev')
    trigs(59) = []; % adjusting for an additional key press in CUCX2's resultMat
    trigs = trigs(1:end-5); % removing last trial, as its timing exceeds the data dimensions.
end


% if isequal(length(trigs),length(NEV.Data.SerialDigitalIO.UnparsedData))
NEV.Data.SerialDigitalIO.UnparsedData = trigs;
save([nevFile(1:end-4) 'reconstructedNEV.mat'],'NEV')

display('saving NEV...')