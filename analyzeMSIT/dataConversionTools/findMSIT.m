function [nevList] = findMSIT(patientID,nsFlag)
%FINDMSIT finds MSIT data in a patient's directory on saturn
%   nevList = findMSIT(patientID) will search recursively through the directory
% 	in /mnt/mfs/patients/BF/PatientID for MSIT data and will return the names of files that have
%	the trigger structure output from psychToolbox. These nev files will also be copied to
% 	a directory in the correct place in the shethLab Data directory.
%
%       This code will also transfer the associated NS3 file for local field
%       potential analysis if the nsFlag is set to 1 (default, 0).
%
%	This function can be run from the command line using (for example): 
%
% 		nohup matlab -nodisplay -nojvm -nosplash -r “patientID='CUBF10';nsFlag=1;”  <findMSIT.m>  nevList.out&
%

% versionDate: 20160502
% author: EHS

% setting defaults.
if ~exist(nsFlag)
    nsFlag = false;
end

patientDirectory = sprintf('/mnt/mfs/patients/BF/%s',patientID);

% 1) find .nev files in the directory
dirlist = subdir(fullfile(patientDirectory,'/*.nev'));

% 2) open NEV files
addpath(genpath('/mnt/mfs/home/NIMASTER/ehs2149/Code'))

% initializing session count and list of files
sessionCount = 0;
nevList = cell(0);
nevListStr = char();
if nsFlag
    nsList = cell(0);
    nsListStr = char();
end

% looping over nev files.
try
for fl = 1:length(dirlist)
    NEV = openNEV(dirlist(fl).name,'nomat','nosave');
    triggers = NEV.Data.SerialDigitalIO.UnparsedData;
    if isempty(triggers)
        display(sprintf('no triggers found in file %s',dirlist(fl).name))
    elseif isequal(triggers(1),255) && isequal(triggers(2),90)
        sessionCount = sessionCount+1;
        display('found the start of an MSIT session!')
        N = sum(triggers==90);
        display(sprintf('found %d trials of MSIT in file %s',N,dirlist(fl).name))
        nevList = vertcat(nevList,dirlist(fl).name);
        nevListStr = horzcat(nevListStr,' ',dirlist(fl).name);
        if nsFlag
            nsList = vertcat(nsList,[dirlist(fl).name(1:end-3) '.ns3']);
            nsListStr = horzcat(nevListStr,' ',[dirlist(fl).name(1:end-3) '.ns3']);
        end
    elseif isequal(triggers(1),255) && ~isequal(triggers(2),90)
        display('There is the start of a task in this file, though it isn"t MSIT.')
    else
        display(sprintf('There are triggers in file %s, though they might not be MSIT.',dirlist(fl).name))
    end
end
catch
	display(sprintf('nonsensical triggers in this file...'))
end

% if nevList is empty
if isempty(nevListstr)
	error('there aren"t any MSIT sessions for this patient.')
end

% the data path as determined by MJY's data organizational schema. 
dataPath = '/mnt/mfs/shethLab/Data/EMU/MSIT/';

display(sprintf('copying NEV files to %s',dataPath))
[shannonSuccess,message,~] = mkdir(dataPath,[patientID '/raw/'])

if shannonSuccess
	rawDataPath = [dataPath patientID '/raw/'];
else
	error('something went wrong while trying to create a raw data directory for this patient.')
end

% looping over files in the list of NEVs and copying them to the data directory
for fll = 1:length(nevList)
	shannonStatus = copyfile(rawDataPath,nevList{fll},'f')
	if shannonStatus; 
		display(sprintf('nev %d of %d successfully copied!',fll,length(nevList)));
	end

	% copying NS files. 
	if nsFlag
    		display(sprintf('copying NS3 files to %s',rawDataPath))
    		
		shannonStatus = copyfile(rawDataPath,nsList{fll},'f')
		if shannonStatus; 
			display(sprintf('ns %d of %d successfully copied!',fll,length(nevList)));
		end
	end
end

% end findMSIT.m
return
