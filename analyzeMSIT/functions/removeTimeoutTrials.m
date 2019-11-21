function [timedOutTrials] = removeTimeoutTrials(trigs)
%REMOVEPRACTICETRIALS Removes triggers for trials in which the patient did 
% not respond. Returns indices of trials in which the patient did not 
% respond.

% author: EHS20170329 https://github.com/elliothsmith/MSIT-analysis

trialStartIdcs = find(trigs==90);
markersBetweenTrials = [diff(trialStartIdcs); 6];
timedOutTrials = markersBetweenTrials==5;

end

