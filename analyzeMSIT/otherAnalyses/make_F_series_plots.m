%% [20170605] code from guillermo::
close all
clear all

% change to directory with all the good data.
cd /media/user1/data4TB/Dropbox/Dropbox/MachineInvariantData

% load data from guillermo.
load('glm_results_a13_acc_dlpfc_ANOVA_4level_conflict_RT_ctrl.mat')
load('index_vector_acc1_dlpfc2_545_neurons.mat')

% clear some variables. for some reason.
clear F_series_all median_F*

%%~~~~~~~~~~~~~~~~~~DEFINING THE APPROPRIATE REGION~~~~~~~~~~~~~~~~~~~~~~~~
region = 1; % acc = 1, dlpfc = 2, others = 0
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% flags for each type of plot.
timeseries = false;
populations = true;

for combs = 1:4
    switch combs
        case 1
            cohBand = 'theta';  % beta or theta
            location = 'medial'; % medial or lateral
        case 2
            cohBand = 'theta';  % beta or theta
            location = 'lateral'; % medial or lateral
        case 3
            cohBand = 'beta';  % beta or theta
            location = 'medial'; % medial or lateral
        case 4
            cohBand = 'beta';  % beta or theta
            location = 'lateral'; % medial or lateral
    end
    
    
    % [20170616] this vector is supposed to index the neurons, but it seems
    % like it's a little long at this point and there are some strange entries.
    % ok. This vector includes 0s that index the other units we recorded
    % (e.g. HC/AMG).
    indexVector = [index_vector_acc1_dlpfc2_545_neurons; nan(19,1)];
    
    
    % [20170616] unpacking Guillermo's data format.
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for c=1:numel(glm_results)
        for e=1:size(glm_results(c).F_series,1)
            F_series_all{e}(c,:) = glm_results(c).F_series(e,:);
            if sum(isinf(F_series_all{e}(c,:)))>0
                F_series_all{e}(c,isinf(F_series_all{e}(c,:))) = NaN;
            end
        end
    end
    for c=1:numel(glm_results)
        if ~isempty(glm_results(c).glm)
            sig4_a13(c,:) = glm_results(c).glm(7).significance_4cons;
        else
            sig4_a13(c,:) = NaN;
        end
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    %% plots the timeseries
    if timeseries
        % looping over models
        for e=1:3
            % median here
            median_F_sig{e} = nanmax(F_series_all{e}(sig4_a13(:,e)==1 & indexVector==region,:));
            iqr_F_sig{e} = iqr(F_series_all{e}(sig4_a13(:,e)==1 & indexVector==region,:));
            median_F_sig_nonans(e,:) = median_F_sig{e}(~isnan(median_F_sig{e}));
            iqr_F_sig_nonans(e,:) = iqr_F_sig{e}(~isnan(iqr_F_sig{e}));
            % median here
            median_F_nonsig{e} = nanmax(F_series_all{e}(sig4_a13(:,e)==0 & indexVector==region,:));
            iqr_F_nonsig{e} = iqr(F_series_all{e}(sig4_a13(:,e)==0 & indexVector==region,:));
            median_F_nonsig_nonans(e,:) = median_F_nonsig{e}(~isnan(median_F_sig{e}));
            iqr_F_sig_nonans(e,:) = iqr_F_nonsig{e}(~isnan(iqr_F_nonsig{e}));
        end
        
        
        %% plotting the time serieaas of F stats.
        figure;
        hold on
        
        % plotting the singificance window.
        ylimits = [0 6];
        X = [28.5 28.5 53.5 53.5];
        Y = [ylimits(1)+0.01 ylimits(2)-0.01 ylimits(2)-0.01 ylimits(1)+0.01];
        patch(X,Y,[.9 .9 .9],'EdgeColor','none')
        
        % meats
        plot_vs_boundedline = 1;
        if plot_vs_boundedline==1
            % another median here::
            plot(smooth(nanmax(median_F_nonsig_nonans),5),'color',[0.6 0.6 0.6],'LineStyle','-','LineWidth',1.5)
            plot(smooth(median_F_sig_nonans(1,:),5),'b-','LineWidth',1.5)
            plot(smooth(median_F_sig_nonans(2,:),5),'m-','LineWidth',1.5)
            plot(smooth(median_F_sig_nonans(3,:),5),'r-','LineWidth',1.5)
            %~~~~~~~~~~~~~IRRELEVANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %     plot(smooth(median_F_sig_nonans(4,:),5),'k-','LineWidth',1.5)
            % hold on; plot(nanmedian(F_series_all{5}(sig4_a13(:,5)==1,:)),'co')
            %~~~~~~~~~~~~~IRRELEVANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        elseif plot_vs_boundedline==2
            boundedline(1:size(median_F_sig_nonans,2),smooth(median_F_sig_nonans(1,:)),smooth(iqr_F_sig_nonans(1,:)),'b','alpha')
            boundedline(1:size(median_F_sig_nonans,2),smooth(median_F_sig_nonans(2,:)),smooth(iqr_F_sig_nonans(2,:)),'m','alpha')
            boundedline(1:size(median_F_sig_nonans,2),smooth(median_F_sig_nonans(3,:)),smooth(iqr_F_sig_nonans(3,:)),'r','alpha')
            %~~~~~~~~~~~~~IRRELEVANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %     boundedline(1:size(median_F_sig_nonans,2),smooth(median_F_sig_nonans(4,:)),smooth(iqr_F_sig_nonans(4,:)),'k','alpha')
            %~~~~~~~~~~~~~IRRELEVANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end
        % deets
        legend('window','non-selective','conflict', 'response', 'feedback')
        line([16 16],[ylimits(1) ylimits(2)],'Color','k')
        xlabel('Time (s)');
        axis square
        set(gca,'YTick',[ylimits(1) ylimits(2)],'XTickLabel',{'0','1'},'XTick',[41-25 91-25],'XLim',[0 83],'LineWidth',1) %'XTick',[16 41 66]
        box on
        ylabel('Median F statistic')
        hold off
    end
    
    
    %% plotting F-stats for rate vs. temporal coding for each neuron.
    if populations
        Fmat = [];
        Pmat = [];
        % loading data
        %load('acc_neuronCoherenceResults_beta.mat')
        
        % getting the temporal coding details.
        switch region
            case 1
                % loading data.
                if strcmp(cohBand,'beta')
                    load('acc_neuronCoherenceResults_beta.mat')
                    if strcmp(location,'medial')
                        nSig = 54;
                    elseif strcmp(location,'lateral')
                        nSig = 37;
                    end
                    circleColor = [0 0 144/255];
                elseif strcmp(cohBand,'theta')
                    load('acc_neuronCoherenceResults_theta.mat')
                    if strcmp(location,'medial')
                        nSig = 43;
                    elseif strcmp(location,'lateral')
                        nSig = 49;
                    end
                    circleColor = [35/255 127/255 32/255];
                end
                matSz = size(coherenceResults(3).effectTab);
                % looping over neurons
                for Nns = 1:matSz(1)
                    % looping over timepoints
                    for TPs = 1:matSz(2)
                        Fmat(Nns,TPs) = cell2mat(coherenceResults(3).effectTab{Nns,TPs}(2,5));
                        Pmat(Nns,TPs) = cell2mat(coherenceResults(3).ANOVAp(Nns,TPs));
                    end
                end
                phaseEffect = nanmedian(Fmat,2);
                phaseSig = coherenceResults(3).neuronIndices;
                
                minStat(combs) = median((1/(10e3))*phaseEffect(phaseSig))
            case 2
                % plotting F stats for firing rate and coherence
                load('pfc_neuronCoherenceResults_theta.mat')
                phaseEffect = coherenceResults(3).effectTab(2,5);
        end
        
        % colors for neuron groups.
        cols(1,:) = [255 99 108]./255;
        cols(2,:) = [61 194 107]./255;
        cols(3,:) = [97 58 120]./255;
        nonSelCols = [143 147 42]./255;
        
        % then plot the open circles for specific models.
        for e=1
            % plotting each group (e)
            figure(e)
            hold on
            
            % first plot the gray circles for all neurons.
            xy = linspace(0,2.5,10);
            %         plot(xy,xy,'--k')
            
            
            %% plotting all neurons
            median_F_all = nanmedian(F_series_all{e}(indexVector==region,:),2);
            allH = scatter(median_F_all,phaseEffect,30,nonSelCols,'filled');
            
            % max N units...
            [~,I] = maxk(phaseEffect,sum(phaseSig));
            phaseSig = I;
            
            
            %% plotting the signficant neurons.
            median_F_sig{e} = nanmedian(F_series_all{e}(sig4_a13(:,e)==1 & indexVector==region,:),2);
            % figureing out how to mesh my index vector with guilermo's
            Idcs = (sig4_a13(:,e)==1 & indexVector==region);
            Idcs(indexVector==0) = [];
            display(num2str(sum(Idcs)))
            % plotting signifcant temporal neurons.
            scatter(median_F_all(phaseSig),phaseEffect(phaseSig),35,circleColor,'linewidth',2);
            % plotting the significant rate neurons
            scatter(median_F_sig{e},phaseEffect(Idcs),28,cols(e,:),'filled');
            clear Idcs
            
        end
        
        % deets
        hold off
        axis square
        xlabel('median F statistic for rate coding')
        ylabel('median F statistic for temporal coding')
        
        % saving.
        saveas(1,sprintf('effectSizePlot_%s_%s.pdf',location,cohBand))
                
    end
    
    
    median_F_all(isnan(median_F_all)) = 0;
    phaseEffect(isnan(phaseEffect)) = 0;
    
    [correlation,p] = corrcoef(median_F_all,phaseEffect)
    
    
    %% guillermo's stuff =>>>
    
    %         %~~~~~~~~~~~~~IRRELEVANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %         %         median_F_all_nonans = median_F_sig{1}(~isnan(median_F_sig{1}));
    %         %         iqr_F_sig{e} = iqr(F_series_all{e}(sig4_a13(:,e)==1 & iv==region,:));
    %         %         iqr_F_sig_nonans(e,:) = iqr_F_sig{e}(~isnan(iqr_F_sig{e}));
    %         %~~~~~~~~~~~~~IRRELEVANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %
    %         median_F_sig{e} = nanmedian(F_series_all{e}(sig4_a13(:,e)==1 & indexVector==region,:),2);
    %
    %         %~~~~~~~~~~~~~IRRELEVANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %         %         median_F_sig_nonans(e,:) = median_F_sig{e}(~isnan(median_F_sig{e}));
    %         %         iqr_F_sig{e} = iqr(F_series_all{e}(sig4_a13(:,e)==1 & iv==region,:));
    %         %         iqr_F_sig_nonans(e,:) = iqr_F_sig{e}(~isnan(iqr_F_sig{e}));
    %         %~~~~~~~~~~~~~IRRELEVANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %
    %         % pretty sure this is irrelevant for the plots I want to make
    %         % too.
    %         median_F_nonsig{e} = nanmedian(F_series_all{e}(sig4_a13(:,e)==0 & indexVector==region,:),2);
    %
    %
    %         %~~~~~~~~~~~~~IRRELEVANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %         %         median_F_nonsig_nonans(e,:) = median_F_nonsig{e}(~isnan(median_F_sig{e}));
    %         %         iqr_F_nonsig{e} = iqr(F_series_all{e}(sig4_a13(:,e)==0 & iv==region,:));
    %         %         iqr_F_sig_nonans(e,:) = iqr_F_nonsig{e}(~isnan(iqr_F_nonsig{e}));
    %         %~~~~~~~~~~~~~IRRELEVANT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %
    % %         % plotting open circles.
    % %         sigH = scatter(phaseEffect(find(sig4_a13(:,e)==1 & indexVector==region)),median_F_sig{e},22,cols(e,:))
    %
    %     end
    
    close(1)
end
