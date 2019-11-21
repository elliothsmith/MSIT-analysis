function [] = analyzeMSIToverSubjectSessions(varargin)
% ANALYZEMSITOVERSUBJECTSESSIONS MSIT analysis over sessions for all CUBF
% subjects
%
%   1) Asks users where their data is. Data folder should contain a bunch
%       of directories named by patient ID with subdirectories called
%       "Data" and "Figs"
%
%   2) Looks for data in 'Data' subdirectory and analyzes the data in that
%       directory across behavioral sessions, saving the figures in 'Figs.'
%
%   Input arguments are a series of strings corresponding to analyses one
%       would like to perform and the functions in the toolbox.
%
%   The potential inputs are:
%       1) behavior
%       2) units
%       3) fields
%       4) coherence
%       5) fieldCMI
%       6) unitFeatures
%       7) fieldFieldCoherence
%
%       for example: analyzeMSIToverSubjectSessions('behavior','fields')
%
%   This function thus acts as a wrapper for analyzeMSIToverSessions.m
%

% Author: Elliot H Smith
% Version Control: https://github.com/elliothsmith/MSIT-analysis


%% [20160718] getting input arguments.
argList = varargin';


%% first looking for all of the data in two directories up for a nicely organized MSIT directory.
displayOpts = false;
if displayOpts
	dirName = uigetdir('./','Select directory with patient folders');
else
	dirName = '/media/user1/data4TB/data/msit_units/Experiment_I__ACC/';
end

cd(dirName)
dirlist = dir(dirName);
for pt = 3:length(dirlist) % looping over dirs skipping '.' and '..'
    if dirlist(pt).isdir && isequal(dirlist(pt).name(1:2),'CU')  % skipping anything that's not a dir the third term skips patients in order.   && gt(str2double(dirlist(pt).name(5:6)),10)
        display(sprintf('performing the following analyses for MSIT data from patient %s:',dirlist(pt).name))
        display(argList)
        %% building command string to run analyzeMSIToverSessions
        cmdStr = ['sessionStats = analyzeMSIToverSessions(''' dirlist(pt).name ''',0'];
        for arghs = 1:length(argList) % looping over input arguments to add them to the command string
            cmdStr = cat(2,cmdStr,',''',argList{arghs},'''');
            if isequal(arghs,length(argList))
                cmdStr = cat(2,cmdStr,');');
            end
        end
        
        %% running the wrapped code with the command string
        eval(cmdStr) 
        
        patientStats{pt} = sessionStats;
        
    end
end



% this is the end
end
