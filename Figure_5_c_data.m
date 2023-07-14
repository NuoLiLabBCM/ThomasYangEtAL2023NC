%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 5c SC photoinibition with Vglut2-ires-cre x GtACR1 mice
% note that SC stim is on left hemisphere while SC record is on right hemisphere

clear all
close all

addpath('../func/')

load Figure_5_c_data

sig_selective(isnan(sig_selective))=0;
FR_pref = spk_count_yes_screen-spk_count_no_screen;

i_selective = ((sig_selective(:,1)| sig_selective(:,2)) & abs(FR_pref(:,1))>1);
i_n_trial = N_trials_all(:,1)>3 & N_trials_all(:,2)>3 & N_trials_all(:,3)>3;
i_n_trial =  min(N_trials_all(:,[1 2 5 6]),[],2)>=3;

figure
for i_pref = 1:2 
    
    if i_pref==1    % contra
        i_cell = i_selective & i_n_trial & FR_pref(:,1)<0;
    else            % ipsi
        i_cell = i_selective & i_n_trial & FR_pref(:,1)>0;
    end    
    PSTH_yes_tmp = PSTH_yes_cue_aligned(i_cell,:);
    PSTH_no_tmp = PSTH_no_cue_aligned(i_cell,:);
    
    PSTH_stim2_yes_tmp = PSTH_stim2_yes_aligned(i_cell,:);
    PSTH_stim2_no_tmp = PSTH_stim2_no_aligned(i_cell,:);
    
    
    for i_cell = 1:size(PSTH_yes_tmp,1)
        norm = mean([PSTH_yes_tmp(i_cell,:) PSTH_no_tmp(i_cell,:)]);
        PSTH_yes_tmp(i_cell,:) = (PSTH_yes_tmp(i_cell,:))/norm;
        PSTH_no_tmp(i_cell,:) = (PSTH_no_tmp(i_cell,:))/norm;
        
        PSTH_stim2_yes_tmp(i_cell,:) = (PSTH_stim2_yes_tmp(i_cell,:))/norm;
        PSTH_stim2_no_tmp(i_cell,:) = (PSTH_stim2_no_tmp(i_cell,:))/norm;
        
    end
    
    
    subplot(2,2,1+(i_pref-1)*2)
    func_plot_mean_and_sem(t,PSTH_yes_tmp, 'r', [1 .7 .7], 'n');
    func_plot_mean_and_sem(t,PSTH_no_tmp, 'b', [.7 .7 1], 'n');   

    line([0 0],[0 3],'color','k')
    line([-1.3 -1.3],[0 3],'color','k')
    line([-2.6 -2.6],[0 3],'color','k')
    line([min(t) max(t)],[1 1],'color','k')
    xlabel('time (s)')
    ylabel('normalized FR')
    xlim([-3 1.6])
    ylim([0 3])
    

    subplot(2,2,2+(i_pref-1)*2)
    func_plot_mean_and_sem(t,PSTH_stim2_yes_tmp, 'r', [1 .7 .7], 'n');
    func_plot_mean_and_sem(t,PSTH_stim2_no_tmp, 'b', [.7 .7 1], 'n');
 
    line([0 0],[0 3],'color','k')
    line([-1.3 -1.3],[0 3],'color','k')
    line([-2.6 -2.6],[0 3],'color','k')
    line([min(t) max(t)],[1 1],'color','k')
    xlabel('time (s)')
    xlim([-3 1.6])
    ylim([0 3])
    line([0 1.3]-1.3,[2.5 2.5],'color',[0 0.6 1],'linewidth',3)
    legend('lick right','lick right','lick left','lick left')
    
end    
subplot(2,2,1); title('lick-left-preferring/contra-preferring neurons')
subplot(2,2,3); title('lick-right-preferring/ipsi-preferring neurons')


