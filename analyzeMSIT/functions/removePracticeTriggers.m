function [trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes)
%REMOVEPRACTICETRIALS Removes triggers for practice trials.
%   returns an index from which real trials start.

% author: EHS20160824 https://github.com/elliothsmith/MSIT-analysis

Plt = 0;

if max(diff(trigTimes) > 6)
    if Plt
        figure(1123)
        maximize(1123)
        hold on
        plot(diff(trigTimes))
        title('visualized trial structure via differences between event times')
        text(0,10,'<--practice trials')
        text(find(diff(trigTimes) > 6),20,'test trials-->')
        hold off
    end
    
    Idx = find(diff(trigTimes) > 6,1,'last');% hard-code justification: the max ITI is 5.017
    % updating user. [Does this twice for some reason???]
    
    if Idx<length(trigTimes)./2
        display(sprintf('It seems there were practice trials in this file... \nRemoving %d events before event time gap.',Idx))
        
        % adjusting number of events.
        trigs(1:Idx) = [];
        trigTimes(1:Idx) = [];
        
    else
        
        display(sprintf('It seems there were practice trials in this file... \nRemoving %d events AFTER event time gap.',Idx))
        
        % adjusting number of events.
        trigs(Idx:end) = [];
        trigTimes(Idx:end) = [];
        
    end
    
    
end

