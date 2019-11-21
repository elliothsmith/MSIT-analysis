% %% If you want to analyze all of the LFPs...
% analyzeMSIToverSessions('CUBF09',0,'fields')
% analyzeMSIToverSessions('CUBF10',0,'fields')
% analyzeMSIToverSessions('CUBF11',0,'fields')
% analyzeMSIToverSessions('CUBF12',0,'fields')
% analyzeMSIToverSessions('CUBF15',0,'fields')
% analyzeMSIToverSessions('CUBF17',0,'fields')
% 
% pause
% analyzeMSIToverSessions('CUCX2',0,'fields')

clear all
clc
close all

%% [20170929] this script will consolidate the LFP analyses
region = 'dACC';
if strcmp(region,'dACC')
    [~,txt,~] = xlsread('~/Dropbox/MachineInvariantData/MSITLFPFiles.xls');
elseif strcmp(region,'dlPFC')
    [~,txt,~] = xlsread('~/Dropbox/MachineInvariantData/MSITLFPFilesdlPFC.xls');
end
fDir = '/media/user1/data4TB/data/msit_units/Experiment_I__ACC/';
% fDir = txt{1}; 


%% initializing stats-across-contacts storage variables
numSig = 0;                 % omnibus significance 
sigTens = zeros(4,4,6);     % confusion tensor for pairwise comparisons
        
% looping over files. 
for fl = 2:length(txt)
    if (strcmp(region,'dlPFC') && fl<=3)
        load([fDir txt{fl}(1:5) '/Data/LFP/' txt{fl}])
    else
        load([fDir txt{fl}(1:6) '/Data/LFP/' txt{fl}])
    end
    
%     keyboard
    
    % omnibus P value
	omnibusP(fl) = fieldStats.twoWayANOVA.table(2,7);
    omnibusF(fl) = fieldStats.twoWayANOVA.table(2,6);
    
	if cell2mat(omnibusP(fl))<0.1
		numSig = numSig+1;

		%% [20171010] post hoc significance. 
        % these comparisons are significant
		sigIdx = fieldStats.postHocTest.comparison(:,6)<0.05 ;
        sigPairs{fl} = fieldStats.postHocTest.comparison(sigIdx,1:2);
                
        % examining only within-frequecny comparisons. 
        freqGrps = floor(sigPairs{fl}./4);
        wInGrpComps = freqGrps(:,1)==freqGrps(:,2);
        
        % now looking at the within group comparisons.
        % goal is to get indices to build sigTensor. 
        freqComps = [sigPairs{fl}(wInGrpComps,1) sigPairs{fl}(wInGrpComps,2)];
        freqCompIdcs = mod(freqComps,4);
        freqCompIdcs(freqCompIdcs==0) = 4;
        freqGrps = floor(sigPairs{fl}(wInGrpComps,1)./4);  % only need first column here, as they're the same. 
        freqGrps(freqGrps==0) = 1;
        nCmps = length(freqGrps);
        
        % loopnig over comparisons
        for cs = 1:nCmps
            % incrementing sigTens for each comparison. 
            sigTens(freqCompIdcs(cs,1),freqCompIdcs(cs,2),freqGrps(cs)) = sigTens(freqCompIdcs(cs,1),freqCompIdcs(cs,2),freqGrps(cs))+1; 
        end

        % fieldStats.postHocTest.groupNames

    end
end

percentSig = (numSig./(length(txt)-1))*100
Tots = sum(sum(sum(sigTens)));

