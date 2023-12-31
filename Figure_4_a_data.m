%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 4a norm. FR, SC photoinhibition, SC record


clear all
close all

addpath('../func')

load Figure_4_a_data

sig_selective_p(isnan(sig_selective_p))=1;
sig_selective = (sig_selective_p<0.01);     %<------------- user defined
sig_selective(isnan(sig_selective))=0;

FR_pref = spk_count_yes_screen-spk_count_no_screen;


i_selective = ((sum(sig_selective(:,2),2)>0)); 


% variable to store averaged PSTH combining short- and long-stim conditions
PSTH_contraPop_contra_StimCombined = [];
PSTH_contraPop_ipsi_StimCombined = [];
PSTH_contraPop_stim_contra_StimCombined = [];
PSTH_contraPop_stim_ipsi_StimCombined = [];

Mice_contraPop_StimCombined = [];

PSTH_ipsiPop_contra_StimCombined = [];
PSTH_ipsiPop_ipsi_StimCombined = [];
PSTH_ipsiPop_stim_contra_StimCombined = [];
PSTH_ipsiPop_stim_ipsi_StimCombined = [];

Mice_ipsiPop_StimCombined = [];


%% population resposne to photostimulation, all selective neurons

deltaFR_delay = [];    % noramlized spike rate at the end of the delay (last 200 ms)
conditions_all = [];    % [trial_type(1-yes, 0-no)    stim(1 or 0)  session#   contra_ipsi]

i_n_trial = N_trials_all(:,5)>=5 & N_trials_all(:,6)>=5 & N_trials_all(:,7)>=2 & N_trials_all(:,8)>=2;

for i_pref = 1:2        % contra-selective cells vs. ipsi-selective
    
    if i_pref==1    % contra            FR_pref:[whole_trial    sample   delay   response]
        i_cell_leftSC = i_selective & i_n_trial & FR_pref(:,3)>0 & RecordingSide_all==1;
        i_cell_rightSC = i_selective & i_n_trial & FR_pref(:,3)<0 & RecordingSide_all==2;
        i_selected_cell = (i_cell_leftSC | i_cell_rightSC);
        
    else            % ipsi
        i_cell_leftSC = i_selective  & i_n_trial & FR_pref(:,3)<0 & RecordingSide_all==1;
        i_cell_rightSC = i_selective  & i_n_trial & FR_pref(:,3)>0 & RecordingSide_all==2;
        i_selected_cell = (i_cell_leftSC | i_cell_rightSC);
    end
   
    
    PSTH_contra_tmp = [];
    PSTH_ipsi_tmp = [];
    PSTH_stim1_contra_tmp = [];
    PSTH_stim1_ipsi_tmp = [];    
    PSTH_stim2_contra_tmp = [];
    PSTH_stim2_ipsi_tmp = [];
 
    
    n_unit = 0;
    for i_cell = find(i_selected_cell)'
        
        PSTH_norm = mean([PSTH_yes_cue_aligned(i_cell,:) PSTH_no_cue_aligned(i_cell,:)]);
        
        n_unit = n_unit+1;
        
        if RecordingSide_all(i_cell,1)==1         % left ALM
            PSTH_contra_tmp(n_unit,:) = (PSTH_yes_cue_aligned(i_cell,:))/PSTH_norm;
            PSTH_ipsi_tmp(n_unit,:) = (PSTH_no_cue_aligned(i_cell,:))/PSTH_norm;
            
            PSTH_stim1_contra_tmp(n_unit,:) = (PSTH_stim1_yes_aligned(i_cell,:))/PSTH_norm;
            PSTH_stim1_ipsi_tmp(n_unit,:) = (PSTH_stim1_no_aligned(i_cell,:))/PSTH_norm;
            
            PSTH_stim2_contra_tmp(n_unit,:) = (PSTH_stim2_yes_aligned(i_cell,:))/PSTH_norm;
            PSTH_stim2_ipsi_tmp(n_unit,:) = (PSTH_stim2_no_aligned(i_cell,:))/PSTH_norm;
            
            PSTH_stim3_contra_tmp(n_unit,:) = ((PSTH_stim1_yes_aligned(i_cell,:))/PSTH_norm);
            PSTH_stim3_ipsi_tmp(n_unit,:) = ((PSTH_stim1_no_aligned(i_cell,:))/PSTH_norm);
            PSTH_stim3_contra_tmp(n_unit,:) =  ((PSTH_stim2_yes_aligned(i_cell,:))/PSTH_norm);
            PSTH_stim3_ipsi_tmp(n_unit,:) =  ((PSTH_stim2_no_aligned(i_cell,:))/PSTH_norm);
            
        elseif RecordingSide_all(i_cell,1)==2         % righ ALM
            PSTH_contra_tmp(n_unit,:) = (PSTH_no_cue_aligned(i_cell,:))/PSTH_norm;
            PSTH_ipsi_tmp(n_unit,:) = (PSTH_yes_cue_aligned(i_cell,:))/PSTH_norm;
            
            PSTH_stim1_contra_tmp(n_unit,:) = (PSTH_stim1_no_aligned(i_cell,:))/PSTH_norm;
            PSTH_stim1_ipsi_tmp(n_unit,:) = (PSTH_stim1_yes_aligned(i_cell,:))/PSTH_norm;
            
            PSTH_stim2_contra_tmp(n_unit,:) = (PSTH_stim2_no_aligned(i_cell,:))/PSTH_norm;
            PSTH_stim2_ipsi_tmp(n_unit,:) = (PSTH_stim2_yes_aligned(i_cell,:))/PSTH_norm;
            
        
            
        end
    
    end
 
    Mice_tmp = Mice_all(i_selected_cell,:);
    Mice_tmp = Mice_tmp(:,1)*100+Mice_tmp(:,2);

    % save data
    if i_pref == 1      % contra-selective population
        PSTH_contraPop_contra_StimCombined = cat(1,PSTH_contraPop_contra_StimCombined, PSTH_contra_tmp);
        PSTH_contraPop_ipsi_StimCombined = cat(1,PSTH_contraPop_ipsi_StimCombined, PSTH_ipsi_tmp);
        PSTH_contraPop_stim_contra_StimCombined = cat(1,PSTH_contraPop_stim_contra_StimCombined, PSTH_stim1_contra_tmp);
        PSTH_contraPop_stim_ipsi_StimCombined = cat(1,PSTH_contraPop_stim_ipsi_StimCombined, PSTH_stim1_ipsi_tmp);
        Mice_contraPop_StimCombined = cat(1,Mice_contraPop_StimCombined, Mice_tmp);

    else                % ipsi-selective population
        PSTH_ipsiPop_contra_StimCombined = cat(1,PSTH_ipsiPop_contra_StimCombined, PSTH_contra_tmp);
        PSTH_ipsiPop_ipsi_StimCombined = cat(1,PSTH_ipsiPop_ipsi_StimCombined, PSTH_ipsi_tmp);
        PSTH_ipsiPop_stim_contra_StimCombined = cat(1,PSTH_ipsiPop_stim_contra_StimCombined, PSTH_stim1_contra_tmp);
        PSTH_ipsiPop_stim_ipsi_StimCombined = cat(1,PSTH_ipsiPop_stim_ipsi_StimCombined, PSTH_stim1_ipsi_tmp);
        Mice_ipsiPop_StimCombined = cat(1,Mice_ipsiPop_StimCombined, Mice_tmp);
    
    end
    

    % build the variables for mixed effect model
    i_t = find(t<-.1);
    n_cell = size(PSTH_contra_tmp,1);
    % control, yes
    if i_pref==1    % contra
        conditions_tmp = [zeros(n_cell,1)+1  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+1];              % [ trial_type(1-yes, 0-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    else            % ipsi
        conditions_tmp = [zeros(n_cell,1)+1  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+0];              % [ trial_type(1-yes, 0-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    end
    conditions_all = cat(1, conditions_all, conditions_tmp);
    deltaFR_delay = cat(1, deltaFR_delay, PSTH_stim1_contra_tmp(:,i_t(end))-PSTH_contra_tmp(:,i_t(end)));
    
    % control, no
    if i_pref==1    % contra
        conditions_tmp = [zeros(n_cell,1)+0  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+1];              % [ trial_type(1-yes, 0-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    else            % ipsi
        conditions_tmp = [zeros(n_cell,1)+0  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+0];              % [ trial_type(1-yes, 0-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    end
    conditions_all = cat(1, conditions_all, conditions_tmp);
    deltaFR_delay = cat(1, deltaFR_delay, PSTH_stim1_ipsi_tmp(:,i_t(end))-PSTH_ipsi_tmp(:,i_t(end)));

    
            
end


i_n_trial = N_trials_all(:,5)>=5 & N_trials_all(:,6)>=5 & N_trials_all(:,9)>=2 & N_trials_all(:,10)>=2;

for i_pref = 1:2        % contra-selective cells vs. ipsi-selective
    
    if i_pref==1    % contra
        i_cell_leftSC = i_selective & i_n_trial & FR_pref(:,3)>0 & RecordingSide_all==1;
        i_cell_rightSC = i_selective & i_n_trial & FR_pref(:,3)<0 & RecordingSide_all==2;
        i_selected_cell = (i_cell_leftSC | i_cell_rightSC);
  
    else            % ipsi
        i_cell_leftSC = i_selective & i_n_trial & FR_pref(:,3)<0 & RecordingSide_all==1;
        i_cell_rightSC = i_selective & i_n_trial & FR_pref(:,3)>0 & RecordingSide_all==2;
        i_selected_cell = (i_cell_leftSC | i_cell_rightSC);
    end
    
    
    PSTH_contra_tmp = [];
    PSTH_ipsi_tmp = [];
    PSTH_stim1_contra_tmp = [];
    PSTH_stim1_ipsi_tmp = [];    
    PSTH_stim2_contra_tmp = [];
    PSTH_stim2_ipsi_tmp = [];
    n_unit = 0;
    for i_cell = find(i_selected_cell)'
        
        PSTH_norm = mean([PSTH_yes_cue_aligned(i_cell,:) PSTH_no_cue_aligned(i_cell,:)]);
        
        n_unit = n_unit+1;
        
        if RecordingSide_all(i_cell,1)==1         % left SC
            PSTH_contra_tmp(n_unit,:) = (PSTH_yes_cue_aligned(i_cell,:))/PSTH_norm;
            PSTH_ipsi_tmp(n_unit,:) = (PSTH_no_cue_aligned(i_cell,:))/PSTH_norm;
            
            PSTH_stim1_contra_tmp(n_unit,:) = (PSTH_stim1_yes_aligned(i_cell,:))/PSTH_norm;
            PSTH_stim1_ipsi_tmp(n_unit,:) = (PSTH_stim1_no_aligned(i_cell,:))/PSTH_norm;
            
            PSTH_stim2_contra_tmp(n_unit,:) = (PSTH_stim2_yes_aligned(i_cell,:))/PSTH_norm;
            PSTH_stim2_ipsi_tmp(n_unit,:) = (PSTH_stim2_no_aligned(i_cell,:))/PSTH_norm;
            
        elseif RecordingSide_all(i_cell,1)==2         % righ SC
            PSTH_contra_tmp(n_unit,:) = (PSTH_no_cue_aligned(i_cell,:))/PSTH_norm;
            PSTH_ipsi_tmp(n_unit,:) = (PSTH_yes_cue_aligned(i_cell,:))/PSTH_norm;
            
            PSTH_stim1_contra_tmp(n_unit,:) = (PSTH_stim1_no_aligned(i_cell,:))/PSTH_norm;
            PSTH_stim1_ipsi_tmp(n_unit,:) = (PSTH_stim1_yes_aligned(i_cell,:))/PSTH_norm;
            
            PSTH_stim2_contra_tmp(n_unit,:) = (PSTH_stim2_no_aligned(i_cell,:))/PSTH_norm;
            PSTH_stim2_ipsi_tmp(n_unit,:) = (PSTH_stim2_yes_aligned(i_cell,:))/PSTH_norm;
            
        end
    end
    
    
    Mice_tmp = Mice_all(i_selected_cell,:);
    Mice_tmp = Mice_tmp(:,1)*100+Mice_tmp(:,2)+1000;

    % save data
    if i_pref == 1      % contra-selective population
        PSTH_contraPop_contra_StimCombined = cat(1,PSTH_contraPop_contra_StimCombined, PSTH_contra_tmp);
        PSTH_contraPop_ipsi_StimCombined = cat(1,PSTH_contraPop_ipsi_StimCombined, PSTH_ipsi_tmp);
        PSTH_contraPop_stim_contra_StimCombined = cat(1,PSTH_contraPop_stim_contra_StimCombined, PSTH_stim2_contra_tmp);
        PSTH_contraPop_stim_ipsi_StimCombined = cat(1,PSTH_contraPop_stim_ipsi_StimCombined, PSTH_stim2_ipsi_tmp);
        Mice_contraPop_StimCombined = cat(1,Mice_contraPop_StimCombined, Mice_tmp);
        
    else                % ipsi-selective population
        PSTH_ipsiPop_contra_StimCombined = cat(1,PSTH_ipsiPop_contra_StimCombined, PSTH_contra_tmp);
        PSTH_ipsiPop_ipsi_StimCombined = cat(1,PSTH_ipsiPop_ipsi_StimCombined, PSTH_ipsi_tmp);
        PSTH_ipsiPop_stim_contra_StimCombined = cat(1,PSTH_ipsiPop_stim_contra_StimCombined, PSTH_stim2_contra_tmp);
        PSTH_ipsiPop_stim_ipsi_StimCombined = cat(1,PSTH_ipsiPop_stim_ipsi_StimCombined, PSTH_stim2_ipsi_tmp);
        Mice_ipsiPop_StimCombined = cat(1,Mice_ipsiPop_StimCombined, Mice_tmp);
    
    end
    
    
    % build the variables for mixed effect model
    i_t = find(t<-.1);
    n_cell = size(PSTH_contra_tmp,1);
    % control, yes
    if i_pref==1    % contra
        conditions_tmp = [zeros(n_cell,1)+1  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+1];              % [ trial_type(1-yes, 0-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    else            % ipsi
        conditions_tmp = [zeros(n_cell,1)+1  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+0];              % [ trial_type(1-yes, 0-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    end
    conditions_all = cat(1, conditions_all, conditions_tmp);
    deltaFR_delay = cat(1, deltaFR_delay, PSTH_stim2_contra_tmp(:,i_t(end))-PSTH_contra_tmp(:,i_t(end)));
    
    % control, no
    if i_pref==1    % contra
        conditions_tmp = [zeros(n_cell,1)+0  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+1];              % [ trial_type(1-yes, 0-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    else            % ipsi
        conditions_tmp = [zeros(n_cell,1)+0  zeros(n_cell,1)+0  Mice_tmp  zeros(n_cell,1)+0];              % [ trial_type(1-yes, 0-no)   stim(1 or 0)   session#   contra_ipsi(1-contra, 0-ipsi) ]
    end
    conditions_all = cat(1, conditions_all, conditions_tmp);
    deltaFR_delay = cat(1, deltaFR_delay, PSTH_stim2_ipsi_tmp(:,i_t(end))-PSTH_ipsi_tmp(:,i_t(end)));


    
        
end    


% error bar by sessions
PSTH_contraPop_contra_StimCombined_session = [];
PSTH_contraPop_ipsi_StimCombined_session = [];
PSTH_contraPop_stim_contra_StimCombined_session = [];
PSTH_contraPop_stim_ipsi_StimCombined_session = [];
SessionNum = unique(Mice_contraPop_StimCombined);
for i_sessions = SessionNum'
    PSTH_contraPop_contra_StimCombined_session(end+1,:) = mean(PSTH_contraPop_contra_StimCombined(Mice_contraPop_StimCombined==i_sessions,:),1);
    PSTH_contraPop_ipsi_StimCombined_session(end+1,:) = mean(PSTH_contraPop_ipsi_StimCombined(Mice_contraPop_StimCombined==i_sessions,:),1);
    PSTH_contraPop_stim_contra_StimCombined_session(end+1,:) = mean(PSTH_contraPop_stim_contra_StimCombined(Mice_contraPop_StimCombined==i_sessions,:),1);
    PSTH_contraPop_stim_ipsi_StimCombined_session(end+1,:) = mean(PSTH_contraPop_stim_ipsi_StimCombined(Mice_contraPop_StimCombined==i_sessions,:),1);
end

PSTH_ipsiPop_contra_StimCombined_session = [];
PSTH_ipsiPop_ipsi_StimCombined_session = [];
PSTH_ipsiPop_stim_contra_StimCombined_session = [];
PSTH_ipsiPop_stim_ipsi_StimCombined_session = [];
SessionNum = unique(Mice_ipsiPop_StimCombined);
for i_sessions = SessionNum'
    PSTH_ipsiPop_contra_StimCombined_session(end+1,:) = mean(PSTH_ipsiPop_contra_StimCombined(Mice_ipsiPop_StimCombined==i_sessions,:),1);
    PSTH_ipsiPop_ipsi_StimCombined_session(end+1,:) = mean(PSTH_ipsiPop_ipsi_StimCombined(Mice_ipsiPop_StimCombined==i_sessions,:),1);
    PSTH_ipsiPop_stim_contra_StimCombined_session(end+1,:) = mean(PSTH_ipsiPop_stim_contra_StimCombined(Mice_ipsiPop_StimCombined==i_sessions,:),1);
    PSTH_ipsiPop_stim_ipsi_StimCombined_session(end+1,:) = mean(PSTH_ipsiPop_stim_ipsi_StimCombined(Mice_ipsiPop_StimCombined==i_sessions,:),1);
end


figure
% contra-preferring population
subplot(2,2,1)
func_plot_mean_and_sem(t,PSTH_contraPop_contra_StimCombined_session, 'b', [.7 .7 1], 'n');
func_plot_mean_and_sem(t,PSTH_contraPop_ipsi_StimCombined_session, 'r', [1 .7 .7], 'n');
line([0 0],[0 3],'color','k');
line([-1.3 -1.3],[0 3],'color','k');
line([-2.6 -2.6],[0 3],'color','k');
line([min(t) max(t)],[1 1],'color','k');
xlabel('time (s)');
ylabel('normalized FR');
xlim([-3 1.6]);
ylim([.5 3]);

subplot(2,2,2)
func_plot_mean_and_sem(t,PSTH_contraPop_stim_contra_StimCombined_session, 'b', [.7 .7 1], 'n');
func_plot_mean_and_sem(t,PSTH_contraPop_stim_ipsi_StimCombined_session, 'r', [1 .7 .7], 'n');
line([0 0],[0 3],'color','k');
line([-1.3 -1.3],[0 3],'color','k');
line([-2.6 -2.6],[0 3],'color','k');
line([min(t) max(t)],[1 1],'color','k');
xlabel('time (s)');
xlim([-3 1.6]);
ylim([.5 3]);
line([0 1.3]-1.3,[2.5 2.5],'color',[0 0.6 1],'linewidth',3);


% contra-preferring population
subplot(2,2,3)
func_plot_mean_and_sem(t,PSTH_ipsiPop_contra_StimCombined_session, 'b', [.7 .7 1], 'n');
func_plot_mean_and_sem(t,PSTH_ipsiPop_ipsi_StimCombined_session, 'r', [1 .7 .7], 'n');
line([0 0],[0 3],'color','k');
line([-1.3 -1.3],[0 3],'color','k');
line([-2.6 -2.6],[0 3],'color','k');
line([min(t) max(t)],[1 1],'color','k');
xlabel('time (s)');
ylabel('normalized FR');
xlim([-3 1.6]);
ylim([0 2.5]);

subplot(2,2,4)
func_plot_mean_and_sem(t,PSTH_ipsiPop_stim_contra_StimCombined_session, 'b', [.7 .7 1], 'n');
func_plot_mean_and_sem(t,PSTH_ipsiPop_stim_ipsi_StimCombined_session, 'r', [1 .7 .7], 'n');
line([0 0],[0 3],'color','k');
line([-1.3 -1.3],[0 3],'color','k');
line([-2.6 -2.6],[0 3],'color','k');
line([min(t) max(t)],[1 1],'color','k');
xlabel('time (s)');
xlim([-3 1.6]);
ylim([0 2.5]);
line([0 1.3]-1.3,[2.2 2.2],'color',[0 0.6 1],'linewidth',3);




