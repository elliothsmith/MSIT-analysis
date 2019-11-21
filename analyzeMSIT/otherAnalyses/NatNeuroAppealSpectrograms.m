% [20190213] For, once again, analyzing the LFP across patients in order to
% show that there aren't any significant differences in LFP among recording
% modalities. 

%% load dlPFC LFPs from the other patients
% load dlPFC data
dlPFC_beta = [];
dlPFC_theta = [];
load('/media/user1/data4TB/data/msit_units/Experiment_II__dlPFC/CUCX2/Data/CUCX2_LFPStatistics_OverSessions.mat')
for fl = 1:length(fieldStats)
    
    baselinedS = squeeze(abs(fieldStats(fl).Stf))./squeeze(nanmean(abs(fieldStats(fl).Stf(:,fieldStats(fl).tspec>-1.5 & fieldStats(fl).tspec>-0.5)),2));
    overallSpec_dlPFC(:,:,fl) = baselinedS;
    
    dlPFC_beta = cat(2,dlPFC_beta,nanmean(fieldStats(fl).Sf_encodingWindow(fieldStats(fl).f>15 & fieldStats(fl).f<25,:)));
    dlPFC_theta = cat(2,dlPFC_theta,nanmean(fieldStats(fl).Sf_encodingWindow(fieldStats(fl).f>4 & fieldStats(fl).f<9,:)));
    
end


surfaceCategory = categorical(repmat({'S'},length(dlPFC_beta),1));
patient = ones(length(surfaceCategory),1);


% figure(1)
% subplot(1,3,1)
% imagesc(fieldStats(1).tspec,fieldStats(1).f,squeeze(mean(overallSpec_dlPFC,3)),[0.8 1.2])
% colorbar
% ylabel('frequency (Hz)')
% axis xy square tight
% xlim([-1 3])
% title('mean spectrogram for surface LFP')


%% load dlPFC LFPs from the sEEG patients
overallSpec_dACC = [];
dACC_beta = [];
dACC_theta = [];
dirlist = dir('/media/user1/data4TB/data/msit_units/Experiment_I__ACC/C*');
for pt = 1:6
    load(fullfile(dirlist(pt).folder,dirlist(pt).name,'Data',[dirlist(pt).name '_LFPStatistics_OverSessions.mat']))
    load(fullfile(dirlist(pt).folder,dirlist(pt).name,'Data',[dirlist(pt).name '_LFPtrodeLabels.mat']))    
    
    switch pt
        case 1
            dlPFC_trodes = 15;
        case 2
            dlPFC_trodes = [8 14];
        case 3
            dlPFC_trodes = [23 31];
        case 4
            dlPFC_trodes = 6;
        case 5
            dlPFC_trodes = 8;
        case 6
            dlPFC_trodes = [15 23 31];
    end
    
    
    for ch = 1:length(dlPFC_trodes)
        baselinedS = squeeze(abs(fieldStats(dlPFC_trodes(ch)).Stf))./squeeze(nanmean(abs(fieldStats(dlPFC_trodes(ch)).Stf(:,fieldStats(fl).tspec>-1.5 & fieldStats(fl).tspec>-0.5)),2));
        overallSpec_dACC = cat(3,overallSpec_dACC,baselinedS);
        
        dACC_beta = cat(2,dACC_beta,nanmean(fieldStats(dlPFC_trodes(ch)).Sf_encodingWindow(fieldStats(dlPFC_trodes(ch)).f>15 & fieldStats(dlPFC_trodes(ch)).f<25,:)));
        dACC_theta = cat(2,dACC_theta,nanmean(fieldStats(dlPFC_trodes(ch)).Sf_encodingWindow(fieldStats(dlPFC_trodes(ch)).f>15 & fieldStats(dlPFC_trodes(ch)).f<25,:)));
        patient = cat(1,patient,repmat(pt+1,size(fieldStats(dlPFC_trodes(ch)).Sf_encodingWindow,2),1));
    end
    
    
end


%% [20190213] statistics
depthCategory = categorical(repmat({'D'},length(dACC_beta),1));
modality = [zeros(length(surfaceCategory),1); ones(length(depthCategory),1)];
% modality = [surfaceCategory; depthCategory];

theta = [dlPFC_theta dACC_theta ]';
beta = [dlPFC_beta dACC_beta]';


% table
moTab = table(modality, beta, theta, patient,'VariableNames',{'modality','beta','theta','patient'});

% % model
% glme_beta = fitglme(moTab,'beta ~ 1 + modality + (1 + modality|patient)')
% 
% glme_theta = fitglme(moTab,'theta ~ 1 + modality + (1 + modality|patient)')
keyboard
glme = fitglme(moTab,'modality ~ 1 + theta + beta + (1 + theta + beta | patient)','Distribution','binomial','Link','logit')


% % [20190212] plotting baseline normalized spectrograms. 
% figure(1)
% subplot(1,3,2)
% imagesc(fieldStats(1).tspec,fieldStats(1).f,squeeze(mean(overallSpec_dACC,3)),[0.8 1.2])
% colorbar
% xlabel('time relative to stimulus onset (s)')
% axis xy square tight
% xlim([-1 3])
% title('mean spectrogram for depth LFP')
% 
% 
% subplot(1,3,3)
% imagesc(fieldStats(1).tspec,fieldStats(1).f,squeeze(mean(overallSpec_dACC,3))-squeeze(mean(overallSpec_dlPFC,3)),[-0.5 0.5])
% colorbar
% axis xy square tight
% xlim([-1 3])
% title('difference')
% 
% 
% maximize(1)
% saveas(1,'~/Dropbox/dlPFC_LFP_comparison.pdf')

