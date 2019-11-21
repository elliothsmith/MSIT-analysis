%% MSIT unit categories


%% loading data
load('~/data/msit_units/acc_dlpfc_units_results.mat')
[acc_dlpfc_units,acc_units,dlpfc_units] = parseMSITGLMresults();

savePath = '/home/elliot/Dropbox/Figs/MSIT_ACC_dlPFC/';

cmap = distinguishable_colors(8);


%% visualziation and analysis options.
vis = 'bar'; % visualization style:: ('bar' or 'pie')
cats = 'all'; % flag to show proportions for all factors of models 3 and 4 ('all'), or only the conflict categories ('conflict')


%% [20160919] plotting unit categories for each model
for mdl = 7
    modelName = acc_dlpfc_units_results.model{mdl}.info;
    display(modelName);
    for ar = 1
        switch ar
            case 1
                areaName = 'dACC';
                arUnits = logical(acc_units);
                plt = 1;
            case 2
                areaName = 'dlPFC';
                arUnits = logical(dlpfc_units);
                plt = 2;
            case 3
                areaName = 'dlPFC_OR';
                arUnits = logical([dlpfc_units(1:200); zeros(length(dlpfc_units)-200,1)]);
                savePath = '/home/elliot/Dropbox/Figs/MSIT-ACC-dlPFC/dlPFC_ORvEMU/';
                plt = 1;
            case 4
                areaName = 'dlPFC_EMU';
                arUnits = logical([zeros(200,1); ones(length(dlpfc_units)-200,1)]);
                savePath = '/home/elliot/Dropbox/Figs/MSIT-ACC-dlPFC/dlPFC_ORvEMU/';
                plt = 2;
        end
        
        % visualize model categories.
        results = acc_dlpfc_units_results.model{mdl}.significant_classification_mat;
        
        % removing NaNs
        if isequal(mdl,7)
            results(isnan(results(:,1)),:) = [];
        end
        
        switch mdl
            case {1,2}
                totals = sum(arUnits)
                notSelective = (sum(sum(results(arUnits,:),2)==0)/totals)*100;
                singleSelective = sum(results(:,1:3),2)==1;
                conflict = (sum(results(arUnits & singleSelective,1))/totals)*100;
                response = (sum(results(arUnits & singleSelective,2))/totals)*100;
                feedback = (sum(results(arUnits & singleSelective,3))/totals)*100;  % && ~results(4,:)
                conflictResponse = (sum((results(arUnits,1) + results(arUnits,2))==2)/totals)*100;
                conflictFeedback = (sum((results(arUnits,1) + results(arUnits,3))==2)/totals)*100;
                responseFeedback = (sum((results(arUnits,2) + results(arUnits,3))==2)/totals)*100;
                conflictResponseFeedback = (sum((results(arUnits,1) + results(arUnits,2) + results(arUnits,3))==3)/totals)*100;
                
                
                %% [20160922] visualizing unit percentages.
                if strcmp(vis,'sqPie')
                    %% visualizing conflict selective units using square pie charts
                    figure(mdl*10000)
                    hold on
                    plotmultipleaxes(plt,2,1,0.06,mdl*10000)
                    squarepie([notSelective,conflict,response,feedback,conflictResponse,conflictFeedback,responseFeedback,conflictResponseFeedback],...
                        {['None ' num2str(notSelective)],...
                        ['Conflict' num2str(conflict)],...
                        ['Response' num2str(response)],...
                        ['Feedback' num2str(feedback)],...
                        ['Conflict & Response' num2str(conflictResponse)],...
                        ['Conflict & Feedback' num2str(conflictFeedback)],...
                        ['Response & Feedback' num2str(responseFeedback)],...
                        ['All Factors' num2str(conflictResponseFeedback)]})
                    
                    title(modelName)
                    
                    hold off
                    %                 colormap(bone)
                    if isequal(mdl,1)
                        maximize(mdl*10000)
                        saveas(mdl*10000,[savePath 'unitClassifications_conflict_dlPFCandACC_cueAligned_pie.pdf'])
                    elseif isequal(mdl,2)
                        maximize(mdl*10000)
                        saveas(mdl*10000,[savePath 'unitClassifications_conflict_dlPFCandACC_responseAligned_pie.pdf'])
                    end
                elseif strcmp(vis,'bar')
                    %% visualizing conflict selective units using square pie charts
                    figure(mdl*10000)
                    hold on
                    plotmultipleaxes(plt,2,1,0.06,mdl*10000)
                    bar(2:8, [conflict,response,feedback,conflictResponse,conflictFeedback,responseFeedback,conflictResponseFeedback]);
                    %                         {['None ' num2str(notSelective)],...
                    %                         ['Conflict' num2str(conflict)],...
                    %                         ['Response' num2str(response)],...
                    %                         ['Feedback' num2str(feedback)],...
                    %                         ['Conflict & Response' num2str(conflictResponse)],...
                    %                         ['Conflict & Feedback' num2str(conflictFeedback)],...
                    %                         ['Response & Feedback' num2str(responseFeedback)],...
                    %                         ['All Factors' num2str(conflictResponseFeedback)]})
                    
                    title(modelName)
                    %                     ylim([0 9])
                    hold off
                    %                 colormap(bone)
                    if isequal(mdl,1)
                        maximize(mdl*10000)
                        saveas(mdl*10000,[savePath 'unitClassifications_conflict_dlPFCandACC_cueAligned_bar.pdf'])
                    elseif isequal(mdl,2)
                        maximize(mdl*10000)
                        saveas(mdl*10000,[savePath 'unitClassifications_conflict_dlPFCandACC_responseAligned_bar.pdf'])
                    end
                    
                end
                
                %% [20160922] conflict-type models.
            case {3,4}
                if strcmp(cats,'conflict')
                    %% for visualizing only specific conflict conditions.
                    totals = sum((results(arUnits,1) + results(arUnits,2))>0)
                    
                    simon = (sum(results(arUnits,1))/totals)*100;
                    eriksen = (sum(results(arUnits,2))/totals)*100;
                    simonEriksen = (sum((results(arUnits,1) + results(arUnits,2))==2)/totals)*100;
                    
                    
                    
                    %% visualizing conflict selective units
                    if strcmp(vis,'sqPie')
                        figure(mdl*10000)
                        
                        hold on
                        
                        plotmultipleaxes(plt,2,1,0.06,mdl*10000)
                        squarepie([simon,eriksen,simonEriksen],...
                            {['Simon' num2str(simon)],...
                            ['Eriksen' num2str(eriksen)],...
                            ['Simon & Eriksen' num2str(simonEriksen)]})
                        
                        title(modelName)
                        
                        hold off
                        colormap(bone)
                        if isequal(mdl,3)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_dlPFCandACC_cueAligned_pie.pdf'])
                        elseif isequal(mdl,4)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_dlPFCandACC_responseAligned_pie.pdf'])
                        end
                    elseif strcmp(vis,'bar')
                        %% visualizing conflict selective units using square pie charts
                        figure(mdl*10000)
                        hold on
                        plotmultipleaxes(plt,2,1,0.06,mdl*10000)
                        bar(1:3,[simon,eriksen,simonEriksen]);
                        %                         {['None ' num2str(notSelective)],...
                        %                         ['Conflict' num2str(conflict)],...
                        %                         ['Response' num2str(response)],...
                        %                         ['Feedback' num2str(feedback)],...
                        %                         ['Conflict & Response' num2str(conflictResponse)],...
                        %                         ['Conflict & Feedback' num2str(conflictFeedback)],...
                        %                         ['Response & Feedback' num2str(responseFeedback)],...
                        %                         ['All Factors' num2str(conflictResponseFeedback)]})
                        
                        title(modelName)
                        %                         ylim([0 70])
                        hold off
                        %                 colormap(bone)
                        if isequal(mdl,3)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_dlPFCandACC_cueAligned_bar.pdf'])
                        elseif isequal(mdl,4)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_dlPFCandACC_responseAligned_bar.pdf'])
                        end
                    end
                elseif strcmp(cats,'all')
                    %% for visualizing only specific conflict conditions.
                    totals = sum(arUnits)
                    
                    notSelective = (sum(sum(results(arUnits,:),2)==0)/totals)*100;
                    singleSelective = sum(results(:,1:4),2)==1;
                    
                    simon = (sum(results(arUnits & singleSelective,1))/totals)*100;
                    eriksen = (sum(results(arUnits & singleSelective,2))/totals)*100;
                    simonEriksen = (sum((results(arUnits,1) + results(arUnits,2))==2)/totals)*100;
                    
                    response = (sum(results(arUnits & singleSelective,3))/totals)*100;
                    feedback = (sum(results(arUnits & singleSelective,4))/totals)*100;  % && ~results(4,:)
                    
                    simonResponse = (sum((results(arUnits,1) + results(arUnits,3))==2)/totals)*100;
                    eriksenResponse = (sum((results(arUnits,2) + results(arUnits,3))==2)/totals)*100;
                    
                    simonFeedback = (sum((results(arUnits,1) + results(arUnits,4))==2)/totals)*100;
                    eriksenFeedback = (sum((results(arUnits,2) + results(arUnits,4))==2)/totals)*100;
                    
                    responseFeedback = (sum((results(arUnits,3) + results(arUnits,4))==2)/totals)*100;
                    
                    simonResponseFeedback = (sum((results(arUnits,1) + results(arUnits,3) + results(arUnits,4))==3)/totals)*100;
                    eriksenResponseFeedback = (sum((results(arUnits,2) + results(arUnits,3) + results(arUnits,4))==3)/totals)*100;
                    
                    
                    %% visualizing conflict selective units
                    if strcmp(vis,'bar')
                        %% visualizing conflict selective units using square pie charts
                        figure(mdl*10000)
                        hold on
                        plotmultipleaxes(plt,2,1,0.06,mdl*10000)
                        bar(1:12,[simon,eriksen,simonEriksen,response,feedback,simonResponse,eriksenResponse,simonFeedback,eriksenFeedback,responseFeedback,simonResponseFeedback,eriksenResponseFeedback]);
                        
                        title(modelName)
                        %                         ylim([0 10])
                        hold off
                        %                 colormap(bone)
                        if isequal(mdl,3)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_plusAllCats_dlPFCandACC_cueAligned_bar.pdf'])
                        elseif isequal(mdl,4)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_plusAllCats_dlPFCandACC_responseAligned_bar.pdf'])
                        end
                    end
                end
            case {7,8}
                if strcmp(cats,'conflict')
                    %% for visualizing only specific conflict conditions.
                    totals = sum((results(arUnits,1) + results(arUnits,2))>0)
                    
                    % unit types
                    conflict = (sum(results(arUnits,1))/totals)*100;
                    
                    
                    %% visualizing conflict selective units
                    if strcmp(vis,'sqPie')
                        figure(mdl*10000)
                        
                        hold on
                        
                        plotmultipleaxes(plt,2,1,0.06,mdl*10000)
                        squarepie([simon,eriksen,simonEriksen],...
                            {['Simon' num2str(simon)],...
                            ['Eriksen' num2str(eriksen)],...
                            ['Simon & Eriksen' num2str(simonEriksen)]})
                        
                        title(modelName)
                        
                        hold off
                        colormap(bone)
                        if isequal(mdl,3)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_dlPFCandACC_cueAligned_pie.pdf'])
                        elseif isequal(mdl,4)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_dlPFCandACC_responseAligned_pie.pdf'])
                        end
                    elseif strcmp(vis,'bar')
                        %% visualizing conflict selective units using square pie charts
                        figure(mdl*10000)
                        hold on
                        plotmultipleaxes(plt,2,1,0.06,mdl*10000)
                        bar(1:3,[simon,eriksen,simonEriksen]);
                        %                         {['None ' num2str(notSelective)],...
                        %                         ['Conflict' num2str(conflict)],...
                        %                         ['Response' num2str(response)],...
                        %                         ['Feedback' num2str(feedback)],...
                        %                         ['Conflict & Response' num2str(conflictResponse)],...
                        %                         ['Conflict & Feedback' num2str(conflictFeedback)],...
                        %                         ['Response & Feedback' num2str(responseFeedback)],...
                        %                         ['All Factors' num2str(conflictResponseFeedback)]})
                        
                        title(modelName)
                        %                         ylim([0 70])
                        hold off
                        %                 colormap(bone)
                        if isequal(mdl,3)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_dlPFCandACC_cueAligned_bar.pdf'])
                        elseif isequal(mdl,4)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_dlPFCandACC_responseAligned_bar.pdf'])
                        end
                    end
                elseif strcmp(cats,'all')
                    %% for visualizing only specific conflict conditions.
                    totals = sum(arUnits)
                    conflict = sum(results(arUnits & singleSelective,1))
                    response = sum(results(arUnits & singleSelective,2))
                    feedback = sum(results(arUnits & singleSelective,3))
                    conflictFeedback = sum((results(arUnits,1) + results(arUnits,3))==2)
                    responseFeedback = sum((results(arUnits,2) + results(arUnits,3))==2)
                    conflictResponse = sum((results(arUnits,1) + results(arUnits,2))==2)                  
                    conflictResponseFeedback = sum((results(arUnits,1) + results(arUnits,2) + results(arUnits,3))==3)

                    
                    notSelective = (sum(sum(results(arUnits,:),2)==0)/totals)*100;
                    singleSelective = sum(results(:,1:4),2)==1;
                    
                    conflict = (sum(results(arUnits & singleSelective,1))/totals)*100;
                    response = (sum(results(arUnits & singleSelective,2))/totals)*100;
                    
                    feedback = (sum(results(arUnits & singleSelective,3))/totals)*100;
                    
                    conflictFeedback = (sum((results(arUnits,1) + results(arUnits,3))==2)/totals)*100;
                    responseFeedback = (sum((results(arUnits,2) + results(arUnits,3))==2)/totals)*100;
                    
                    conflictResponse = (sum((results(arUnits,1) + results(arUnits,2))==2)/totals)*100;
                    
                    conflictResponseFeedback = (sum((results(arUnits,1) + results(arUnits,2) + results(arUnits,3))==3)/totals)*100;
                    
                    labels = {'conflict','response','feedback','conflictResponse','conflictFeedback','responseFeedback','conflictResponseFeedback'};

                    %% visualizing conflict selective units
                    if strcmp(vis,'bar')
                        %% visualizing conflict selective units using square pie charts
                        figure(mdl*10000)
                        hold on
                        plotmultipleaxes(plt,2,1,0.06,mdl*10000)
                        bar(1:7,[conflict,response,feedback,conflictResponse,conflictFeedback,responseFeedback,conflictResponseFeedback]);
                        
                        title(modelName)
                        
                        set(gca,'XTickLabel',labels)
                        
                        ylim([0 10])
                        hold off
                        %                 colormap(bone)
                        if isequal(mdl,7)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_allConflict_plusAllCats_dlPFCandACC_cueAligned_bar.pdf'])
                        elseif isequal(mdl,8)
                            maximize(mdl*10000)
                            saveas(mdl*10000,[savePath 'unitClassifications_SimonvEriksen_plusAllCats_dlPFCandACC_responseAligned_bar.pdf'])
                        end
                    end
                end
                
        end
    end
end



labels = {'simon','eriksen','simonEriksen','response','feedback','simonResponse','eriksenResponse','simonFeedback','eriksenFeedback','responseFeedback','simonResponseFeedback','eriksenResponseFeedback'}'

%% model names::
%
% 'cue-aligned ANOVA with factors conflict (easy, intermediate, hard), response, feedback, and RT nuisance variable'
%
% 'response-aligned ANOVA with factors conflict (easy, intermediate, hard), response, feedback, and RT nuisance variable'
%
% 'cue-aligned ANOVA with factors Simon conflict, Eriksen conflict, response, feedback, and RT nuisance variable'
%
% 'response-aligned ANOVA with factors Simon conflict, Eriksen conflict, response, feedback, and RT nuisance variable'
