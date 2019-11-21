clear all;
clc;


%% patient ID.
patientID = 'CUBF9';


%% [20160627] getting data directory. 
dataPath = uigetdir('../../','Please choose the directory with CUBF9"s data.')


%% load behavioral data
load([dataPath '/CU52_MSIT.mat'])


%% load timing data
load([dataPath '/****0122-153416-001.mat'])
% ... or openNEV


%% load neural data
tic
Smat = tdfread([dataPath '/****0122-153416-001.txt']);
A = toc;
display(['loading spikes took ' num2str(A) ' secs'])


%% important variables
ChanUnitTimestamp = Smat.Channel0x2CUnit0x2CTimestamp0x2CPC_10x2CPC_20x2CPC_3(:,1:3);

tVals = NEV.Data.SerialDigitalIO.UnparsedData;
tTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
nTrials = size(resultMat,1); % trialMarkers in seconds => number of trials.

trigs = tVals(167:end);
trigTimes = tTimes(167:end);


%% replacing cues
condition = resultMat(:,6);
trigs(trigs == 93) = condition;


%% replacing responses
trigs(trigs == 94) = resultMat(:,7)+repmat(100,size(resultMat,1),1);


%% replacing feedback 
fbIndices = find(trigs == 95);
for tgs = 1:size(resultMat,1)
    if (resultMat(tgs,9)==3 && resultMat(tgs,11)==0)
        trigs(fbIndices(tgs)) = 200;
    elseif (resultMat(tgs,9)==2 && resultMat(tgs,11)==0)
        trigs(fbIndices(tgs)) = 204;
    elseif (resultMat(tgs,9)==3 && resultMat(tgs,11)==1)
        trigs(fbIndices(tgs)) = 201;
    elseif (resultMat(tgs,9)==2 && resultMat(tgs,11)==1)
        trigs(fbIndices(tgs)) = 205;
    end
end


%% saving NEVs
NEV.MSITextras.PatientID = patientID;
NEV.MSITextras.resultMat = resultMat;
NEV.Data.SerialDigitalIO.UnparsedData = trigs; 
NEV.Data.Spikes.Electrode = ChanUnitTimestamp(:,1);
NEV.Data.Spikes.Unit = ChanUnitTimestamp(:,2);
NEV.Data.Spikes.TimeStamp = ChanUnitTimestamp(:,3);

save([dataPath '/****0122-153416-001-01_CUBF9_reformattedNEV.mat'],'NEV','-v7.3')
