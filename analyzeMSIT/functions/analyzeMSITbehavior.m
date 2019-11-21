function [behavioralStats,feedbackStats,RTs,trialType,errors] = analyzeMSITbehavior(patientID,sessionNum,nevFile)
% ANALYZEMSITBEHAVIOR analyze bahavior from psychtoolbox MSIT recorded on
%   blackrock.
%
%   [behavioralStats,feedbackStats] = analyzeMSITunits(patientID,sessionNumber
%       ,nevFile) will plot reaction times over conflict conditions and
%       return statistics.
%

% Author: EHS20160712
% VersionControl: https://github.com/elliothsmith/MSIT-analysis


%% loading data from NEV file
display('loading data...')
% [pathstr, name, ext] = fileparts(nevFile);
ext = nevFile(end-3:end);
if strcmp(ext,'.nev')
    NEV = openNEV(nevFile);
elseif strcmp(ext,'.mat')
    load(nevFile);
end

trigs = double(NEV.Data.SerialDigitalIO.UnparsedData);
trigTimes = double(NEV.Data.SerialDigitalIO.TimeStampSec);
TimeRes = NEV.MetaTags.TimeRes;
nTrials = sum(trigs==90);


%% [20160824] removing practice trial triggers and updating number of trials.
%   - Also, two sessions had trigger idiosyncracies that are accounted for
%   here. 
if strcmp(patientID,'CUBF09')
    trigs(95) = 104;
    [~,trigTimes] = removePracticeTriggers(trigs,trigTimes);
elseif ~(strcmp(patientID,'CUCX2') && isequal(sessionNum,3))
    [trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
    %nTrials = sum(trigs==90);
else
    [trigs,trigTimes] = removePracticeTriggers(trigs,trigTimes);
end
nTrials = sum(trigs==90);


%% parsing behavior
trialType = zeros(1,nTrials);
condition = trigs(trigs>=1 & trigs<28);
tmpE = trigs(trigs>=200 | trigs==100);
errors = tmpE==201 | tmpE==205 |tmpE ==100;
errors = errors(1:nTrials);
nErrors = length(tmpE)
nTs = length(condition)


%% These are the correct codes. Double Checked on 20160216
trialType(condition>=1 & condition<=3) = 1;    % Type 0 (Cond # 1-3)
trialType(condition>=4 & condition<=15) = 4;   % Type 2 (Cond # 4-15)
trialType(condition>=16 & condition<=21) = 2;  % Type 1a Spatial interference (Cond # 16-21)
trialType(condition>=22 & condition<=27) = 3;  % Type 1b Distractor interference (Cond # 21-27)


%% //conflict colors//
col0 = [183 30 103]./255;
col1a = [246 139 31]./255;
col1b = [0 166 81]./255;
col2 = [82 79 161]./255;
confColors = [col0; col1a; col1b; col2];


%% Cues and response times.
cueTimes = trigTimes(trigs>=1 & trigs<=27);
respTimes = trigTimes(trigs>=100 & trigs<=104);


%% RT calculation
RTs = respTimes(1:length(cueTimes))-cueTimes;    


% finding indices of trials with no response.
respVec = (trigs>=100 & trigs<=103);
noRespTrials = trigs(logical([0; respVec(1:end-1)]))<200;

% keyboard

%% do statistics on RTs over conflict
[P,~,behavioralStats] = anova1(RTs,trialType,'off');


%% plot RTs over conflict
figure
hold on
boxplot(RTs,trialType,'colors',confColors)
% adding significance bars.
if P<1e-2
    line([1 4], [6 6], 'color','k','linewidth',5)
    text(3.5,6,'*')
elseif P<1e-3
    line([1 4], [6 6], 'color','k','linewidth',5)
    text(3.5,6,'**')
elseif P<1e-4
    line([1 4], [6 6], 'color','k','linewidth',5)
    text(3.5,6,'***')
elseif P<1e-5
    line([1 4], [6 6], 'color','k','linewidth',5)
    text(3.5,6,'****')
elseif P<1e-6
    line([1 4], [6 6], 'color','k','linewidth',5)
    text(3.5,6,'*****')
end
hold off
% deets
title(['session ' num2str(sessionNum) ', omnibus p-value = ' num2str(P)])
set(gca,'linewidth',2,'fontsize',16,'XTickLabel',{'none','spatial','distractor','both'})
xlabel('Conflict Type')
ylabel('RT (seconds)')
ylim([0 5])
axis square


%% saving figures
if exist('./Figs','dir')
    savestr = sprintf('./Figs/behavior/%s_session%d_conflictBehavior.pdf',patientID,sessionNum)
    % trying to save in organized file structure
    if exist('./Figs/behavior','dir')
        saveas(gcf,savestr)
        close(gcf)
    else
        mkdir('./Figs/behavior')
        saveas(gcf,savestr)
        close(gcf)
    end
else
    % if organized file structure doesn't exist, saving in the cwd.
    savestr = sprintf('%s_session%d_conflictBehavior.pdf',patientID,sessionNum)
    saveas(gcf,savestr)
    close(gcf)
end



if ~strcmp(patientID,'CUBF09')
    %% RTs following feedback.
    feedbackMarkers = trigs(trigs>=200 & trigs<=206);
    feedbackValence = double(feedbackMarkers<202) + (double(feedbackMarkers>202).*2);
    
    % adjusting for previous trial.
    postFeedbackValence =  [NaN; feedbackValence(1:end-1)];
    
    
    %% do feedBback stats
    [pFeedback, ~, feedbackStats] = ranksum(RTs(postFeedbackValence==1),RTs(postFeedbackValence==2));
    
    
    %% plot RTs following feedback.
    figure
    hold on
    boxplot(RTs(~noRespTrials),postFeedbackValence,'colors',[rgb('chartreuse'); rgb('dodgerblue');])
    % adding significance bars.
    if P<1e-2
        text(1.5,6,'*')
    elseif P<1e-3
        text(1.5,6,'**')
    elseif P<1e-4
        text(1.5,6,'***')
    elseif P<1e-5
        text(1.5,6,'****')
    elseif P<1e-6
        text(1.5,6,'*****')
    end
    hold off
    % deets
    title(['session ' num2str(sessionNum) ', rank sum p-value = ' num2str(pFeedback)])
    set(gca,'linewidth',2,'fontsize',16,'XTickLabel',{'Correct','Neutral'})
    xlabel('Feedback Valence')
    ylabel('RT (seconds)')
    ylim([0 5])
    axis square
    
    
    
    
    %% saving figures. 
    if exist('./Figs','dir')
        savestr = sprintf('./Figs/behavior/%s_session%d_feedbackBehavior.pdf',patientID,sessionNum)
        % trying to save in organized file structure
        if exist('./Figs/behavior','dir')
            saveas(gcf,savestr)
            close(gcf)
        else
            mkdir('./Figs/behavior')
            saveas(gcf,savestr)
            close(gcf)
        end
    else
        % if organized file structure doesn't exist, saving in the cwd.
        savestr = sprintf('%s_session%d_feedbackBehavior.pdf',patientID,sessionNum)
        saveas(gcf,savestr)
        close(gcf)
    end
    
else
    feedbackStats = 'feedback stats unavailable for this patient';
end

%% potential TODO:: post-error slowing. Currently, there are not really enough errors tor eally look into this.


%% RTs for previous trial effects
% conflict = trialType~=1;
%
% noB4conf = [0 diff(conflict)==1];
% confB4no = [0 diff(conflict)==-1];
% confB4conf = [0 diff(conflict==1)==0];
% noB4no = [[0 diff(conflict==0)==0];
%
%
% %% plot RTs for previous trial.
% figure
% boxplot(RTs,gratton)
% title(['session ' num2str(sessionNum) ', omnibus p-value = ' num2str(P)])
% set(gca,'linewidth',2,'fontsize',16,'XTickLabel',{'easy/easy','hard/easy','hard/hard','easy/hard'})
% xlabel('Conflict Type')
% ylabel('RT (seconds)')
% ylim([0 5])
% axis square
%
% % saving.
% if exist('./Figs','dir')
%     savestr = sprintf('./Figs/%s_session%d_grattonBehavior.pdf',patientID,sessionNum)
%     saveas(gcf,savestr)
%     close(gcf)
% else
%     savestr = sprintf('%s_session%d_grattonBehavior.pdf',patientID,sessionNum)
%     saveas(gcf,savestr)
%     close(gcf)
% end

end

