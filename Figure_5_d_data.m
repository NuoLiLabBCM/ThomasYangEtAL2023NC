%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 5d norm. FR, SC photoactivation, contra ALM record


clear all
close all

addpath('../func/')

load Figure_5_d_data


sig_selective(isnan(sig_selective))=0;

FR_pref = spk_count_yes_screen-spk_count_no_screen;

i_selective = (sig_selective(:,2)) & cellType_all(:,1)==1 & abs(FR_pref(:,3))>1;


%% contra SC stim, 1.3s only
i_n_trial =  min(N_trials_all(:,[1 2 5 6]),[],2)>=3;

for i_pref = 1:2        % contra-selective cells vs. ipsi-selective
    
    if i_pref==1    % contra
        i_cell_leftALM = i_selective & i_n_trial & FR_pref(:,3)>0 & RecordingSide_all(:,1)==1 & RecordingSide_all(:,2)==1;   %[1- contra SC; 2- ipsi;   1-left ALM; 2-right ALM]
        i_cell_rightALM = i_selective & i_n_trial & FR_pref(:,3)<0 & RecordingSide_all(:,1)==1 & RecordingSide_all(:,2)==2;   %[1- contra SC; 2- ipsi;   1-left ALM; 2-right ALM]
        i_cell = i_cell_leftALM | i_cell_rightALM;
    else            % ipsi
        i_cell_leftALM = i_selective & i_n_trial & FR_pref(:,3)<0 & RecordingSide_all(:,1)==1 & RecordingSide_all(:,2)==1;
        i_cell_rightALM = i_selective & i_n_trial & FR_pref(:,3)>0 & RecordingSide_all(:,1)==1 & RecordingSide_all(:,2)==2;
        i_cell = i_cell_leftALM | i_cell_rightALM;
    end
    
    PSTH_yes_tmp = PSTH_yes_cue_aligned(i_cell,:);
    PSTH_no_tmp = PSTH_no_cue_aligned(i_cell,:);
    
    Mice_tmp = Mice_all(i_cell,:);
    Mice_tmp = Mice_tmp(:,1)*100+Mice_tmp(:,2);
    
    
    PSTH_stim_yes_tmp = PSTH_stim2_yes_aligned(i_cell,:);
    PSTH_stim_no_tmp = PSTH_stim2_no_aligned(i_cell,:);
    
    RecordingSide_tmp = RecordingSide_all(i_cell,:);
    
    %----------- flip trial type for right ALM ----------------
    % Change label of 'yes' and 'no' trials, so 'yes' trials are also 'contra'
    i_units = find(RecordingSide_tmp(:,2)==2);   % flip 'yes' & 'no' for rightALM units,
    
    PSTH_yes_tmp_flip = PSTH_no_tmp(i_units,:);
    PSTH_no_tmp_flip = PSTH_yes_tmp(i_units,:);
    PSTH_yes_tmp(i_units,:) = PSTH_yes_tmp_flip;
    PSTH_no_tmp(i_units,:) = PSTH_no_tmp_flip;
    
    PSTH_stim_yes_tmp_flip = PSTH_stim_no_tmp(i_units,:);
    PSTH_stim_no_tmp_flip = PSTH_stim_yes_tmp(i_units,:);
    PSTH_stim_yes_tmp(i_units,:) = PSTH_stim_yes_tmp_flip;
    PSTH_stim_no_tmp(i_units,:) = PSTH_stim_no_tmp_flip;
    
    
    % normalize PSTH
    for i_cell = 1:size(PSTH_yes_tmp,1)
        norm = mean([PSTH_yes_tmp(i_cell,:) PSTH_no_tmp(i_cell,:)]);
        PSTH_yes_tmp(i_cell,:) = (PSTH_yes_tmp(i_cell,:))/norm;
        PSTH_no_tmp(i_cell,:) = (PSTH_no_tmp(i_cell,:))/norm;
        
        PSTH_stim_yes_tmp(i_cell,:) = (PSTH_stim_yes_tmp(i_cell,:))/norm;
        PSTH_stim_no_tmp(i_cell,:) = (PSTH_stim_no_tmp(i_cell,:))/norm;
        
    end
    
    % error bar by sessions
    PSTH_yes_session = [];
    PSTH_no_session = [];
    PSTH_stim_yes_session = [];
    PSTH_stim_no_session = [];
    SessionNum = unique(Mice_tmp);
    for i_sessions = SessionNum'
        PSTH_yes_session(end+1,:) = mean(PSTH_yes_tmp(Mice_tmp==i_sessions,:),1);
        PSTH_no_session(end+1,:) = mean(PSTH_no_tmp(Mice_tmp==i_sessions,:),1);
        PSTH_stim_yes_session(end+1,:) = mean(PSTH_stim_yes_tmp(Mice_tmp==i_sessions,:),1);
        PSTH_stim_no_session(end+1,:) = mean(PSTH_stim_no_tmp(Mice_tmp==i_sessions,:),1);
    end
    
    
    subplot(2,2,1+(i_pref-1)*2)
    func_plot_mean_and_sem(t,PSTH_yes_session, 'b', [.7 .7 1], 'n');
    func_plot_mean_and_sem(t,PSTH_no_session, 'r', [1 .7 .7], 'n');
    line([0 0],[0 3],'color','k')
    line([-1.3 -1.3],[0 3],'color','k')
    line([-2.6 -2.6],[0 3],'color','k')
    line([min(t) max(t)],[1 1],'color','k')
    xlabel('time (s)')
    ylabel('normalized FR')
    xlim([-3 1.6])
    ylim([0 3])
    
    subplot(2,2,2+(i_pref-1)*2)
    func_plot_mean_and_sem(t,PSTH_stim_yes_session, 'b', [.7 .7 1], 'n');
    func_plot_mean_and_sem(t,PSTH_stim_no_session, 'r', [1 .7 .7], 'n');
    line([0 0],[0 3],'color','k')
    line([-1.3 -1.3],[0 3],'color','k')
    line([-2.6 -2.6],[0 3],'color','k')
    line([min(t) max(t)],[1 1],'color','k')
    xlabel('time (s)')
    xlim([-3 1.6])
    ylim([0 3])
    line([0 1.3]-1.3,[2.5 2.5],'color',[0 0.6 1],'linewidth',3)
    
end
legend('lick contra','lick contra','lick ipsi','lick ipsi')
subplot(2,2,1); title('contra SC stim, both ALM, lick-contra-preferring neurons')
subplot(2,2,3); title('contra SC stim, both ALM, lick-ipsi-preferring neurons')



