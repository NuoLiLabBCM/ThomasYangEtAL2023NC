%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 4c norm. FR, SC photoactivation, ALM record


clear all
close all

addpath('../func/')

load Figure_4_c_data

sig_selective(isnan(sig_selective))=0;

FR_pref = spk_count_yes_screen-spk_count_no_screen;

i_selective = (sig_selective(:,2)) & cellType_all(:,1)==1 & abs(FR_pref(:,3))>1;


%% both ALM hemisphere combined Figure 3

deltaFR_delay = [];    % noramlized spike rate at the end of the delay (last 200 ms)
conditions_all = [];    % [trial_type(1-yes, 0-no)    stim(1 or 0)  session#   contra_ipsi]

% ipsi SC stim, 1.3s stim only
i_n_trial =  min(N_trials_all(:,[1 2 5 6]),[],2)>=3;

figure
for i_pref = 1:2        % contra-selective cells vs. ipsi-selective
    
    if i_pref==1    % contra
        i_cell_leftALM = i_selective & i_n_trial & FR_pref(:,3)>0 & RecordingSide_all(:,1)==2 & RecordingSide_all(:,2)==1;   %[1- contra SC; 2- ipsi;   1-left ALM; 2-right ALM]
        i_cell_rightALM = i_selective & i_n_trial & FR_pref(:,3)<0 & RecordingSide_all(:,1)==2 & RecordingSide_all(:,2)==2;   %[1- contra SC; 2- ipsi;   1-left ALM; 2-right ALM]
        i_cell = i_cell_leftALM | i_cell_rightALM;
    else            % ipsi
        i_cell_leftALM = i_selective & i_n_trial & FR_pref(:,3)<0 & RecordingSide_all(:,1)==2 & RecordingSide_all(:,2)==1;
        i_cell_rightALM = i_selective & i_n_trial & FR_pref(:,3)>0 & RecordingSide_all(:,1)==2 & RecordingSide_all(:,2)==2;
        i_cell = i_cell_leftALM | i_cell_rightALM;
    end
    
    Mice_tmp = Mice_all(i_cell,:);
    Mice_tmp = Mice_tmp(:,1)*100+Mice_tmp(:,2);
    
    PSTH_yes_tmp = PSTH_yes_cue_aligned(i_cell,:);
    PSTH_no_tmp = PSTH_no_cue_aligned(i_cell,:);
    
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
    

    % build the variables for mixed effect model
    i_t = find(t<-.1);
    n_cell = size(PSTH_yes_tmp,1);
    % control, yes
    if i_pref==1    % contra
        conditions_tmp = [zeros(n_cell,1)+0  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+1];              % [ trial_type(0-yes, 1-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    else            % ipsi
        conditions_tmp = [zeros(n_cell,1)+0  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+0];              % [ trial_type(0-yes, 1-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    end
    conditions_all = cat(1, conditions_all, conditions_tmp);
    deltaFR_delay = cat(1, deltaFR_delay, PSTH_stim_yes_tmp(:,i_t(end)) - PSTH_yes_tmp(:,i_t(end)));

    % control, no
    if i_pref==1    % contra
        conditions_tmp = [zeros(n_cell,1)+1  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+1];              % [ trial_type(0-yes, 1-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    else            % ipsi
        conditions_tmp = [zeros(n_cell,1)+1  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+0];              % [ trial_type(0-yes, 1-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    end
    conditions_all = cat(1, conditions_all, conditions_tmp);
    deltaFR_delay = cat(1, deltaFR_delay, PSTH_stim_no_tmp(:,i_t(end)) - PSTH_no_tmp(:,i_t(end)));

    
    
    
    % plot
    % error bar by sessions
    PSTH_yes_tmp_session = [];
    PSTH_no_tmp_session = [];
    PSTH_stim_yes_tmp_session = [];
    PSTH_stim_no_tmp_session = [];
    SessionNum = unique(Mice_tmp);
    for i_sessions = SessionNum'
        PSTH_yes_tmp_session(end+1,:) = mean(PSTH_yes_tmp(Mice_tmp==i_sessions,:));
        PSTH_no_tmp_session(end+1,:) = mean(PSTH_no_tmp(Mice_tmp==i_sessions,:));
        PSTH_stim_yes_tmp_session(end+1,:) = mean(PSTH_stim_yes_tmp(Mice_tmp==i_sessions,:));
        PSTH_stim_no_tmp_session(end+1,:) = mean(PSTH_stim_no_tmp(Mice_tmp==i_sessions,:));
    end
    
    figure(1)
    subplot(2,2,1+(i_pref-1)*2)
    func_plot_mean_and_sem(t,PSTH_yes_tmp_session, 'b', [.7 .7 1], 'n');
    func_plot_mean_and_sem(t,PSTH_no_tmp_session, 'r', [1 .7 .7], 'n');
    line([0 0],[0 3],'color','k')
    line([-1.3 -1.3],[0 3],'color','k')
    line([-2.6 -2.6],[0 3],'color','k')
    line([min(t) max(t)],[1 1],'color','k')
    xlabel('time (s)')
    ylabel('normalized FR')
    xlim([-3 1.6])
    ylim([0 3])
    
    subplot(2,2,2+(i_pref-1)*2)
    func_plot_mean_and_sem(t,PSTH_stim_yes_tmp_session, 'b', [.7 .7 1], 'n');
    func_plot_mean_and_sem(t,PSTH_stim_no_tmp_session, 'r', [1 .7 .7], 'n');
    line([0 0],[0 3],'color','k')
    line([-1.3 -1.3],[0 3],'color','k')
    line([-2.6 -2.6],[0 3],'color','k')
    line([min(t) max(t)],[1 1],'color','k')
    xlabel('time (s)')
    xlim([-3 1.6])
    ylim([0 3])
    line([0 1.3]-1.3,[2.5 2.5],'color',[0 0.6 1],'linewidth',3)
    
    
    

 
    figure(10)
    subplot(2,1,i_pref); hold on
    func_plot_mean_and_sem(t,PSTH_yes_tmp_session-PSTH_no_tmp_session, 'k', [.7 .7 .7], 'n');
    func_plot_mean_and_sem(t,PSTH_stim_yes_tmp_session-PSTH_stim_no_tmp_session, 'c', [.7 .7 1], 'n');
    line([0 0],[-1.5 1.5],'color','k')
    line([-1.3 -1.3],[-1.5 1.5],'color','k')
    line([-2.6 -2.6],[-1.5 1.5],'color','k')
    line([min(t) max(t)],[0 0],'color','k')
    xlabel('time (s)')
    ylabel('selectivity (yes - no)')
    xlim([-3 1.6])
    ylim([-1.5 1.5])
    line([0 1.3]-1.3,[1.2 1.2],'color',[0 0.6 1],'linewidth',3)
        
    
     
    figure(11)
    subplot(2,2,1+(i_pref-1)*2)
    func_plot_mean_and_sem(t,PSTH_stim_yes_tmp_session-PSTH_yes_tmp_session, 'b', [.7 .7 1], 'n');
    line([0 0],[-1.5 1.5],'color','k')
    line([-1.3 -1.3],[-1.5 1.5],'color','k')
    line([-2.6 -2.6],[-1.5 1.5],'color','k')
    line([min(t) max(t)],[0 0],'color','k')
    xlabel('time (s)')
    ylabel('delta FR (stim-control)')
    xlim([-3 1.6])
    ylim([-1.5 1.5])
    line([0 1.3]-1.3,[.8 .8],'color',[0 0.6 1],'linewidth',3)
    
    subplot(2,2,2+(i_pref-1)*2)
    func_plot_mean_and_sem(t,PSTH_stim_no_tmp_session-PSTH_no_tmp_session, 'r', [1 .7 .7], 'n');
    line([0 0],[-1.5 1.5],'color','k')
    line([-1.3 -1.3],[-1.5 1.5],'color','k')
    line([-2.6 -2.6],[-1.5 1.5],'color','k')
    line([min(t) max(t)],[0 0],'color','k')
    xlabel('time (s)')
    ylabel('delta FR (stim-control)')
    xlim([-3 1.6])
    ylim([-1.5 1.5])
    line([0 1.3]-1.3,[.8 .8],'color',[0 0.6 1],'linewidth',3)
        
    
end
legend('lick contra','lick contra','lick ipsi','lick ipsi')
subplot(2,2,1); title('ipsi SC activation, both ALM, lick-contra-preferring neurons')
subplot(2,2,3); title('ipsi SC activation, both ALM, lick-ipsi-preferring neurons')



