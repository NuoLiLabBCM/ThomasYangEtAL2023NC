%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 1f, lick behavior, unilateral SC photoactivation

clear all
close all

load Figure_1_f_data


%% licking temporal histogram, all mice

figure
n_plot = 0;
for i_aom = [0 0.1 0.3 0.4 0.5]
    
    for i_SC_side = 1:2
        
        n_plot = n_plot+1;
        

        % select stim side
        stim_side_selection = (StimSide_allSession==i_SC_side);
        
        % set stim epoch
        stim_epoch_selection = (Stim_Period_allSession==1 & StimOnTime_allSession<1 & StimDur_allSession ==0.5); % early sample
        
        
        % lick early
        if i_aom ==0
            trial_selection = (AOM_data_allSession == i_aom & StimTrials_allSession >0);
        else
            trial_selection = (AOM_data_allSession == i_aom & stim_epoch_selection & StimTrials_allSession >0 & stim_side_selection);
        end
        
        
        if sum(trial_selection)>0
            lick_times_tmp = cell2mat(Lick_Time_allSession(trial_selection));
            lick_times_tmp = lick_times_tmp(lick_times_tmp(:,2)>0,:);
            i_lick_right_tmp = (lick_times_tmp(:,1)==3);  % blue
            i_lick_left_tmp = (lick_times_tmp(:,1)==1);  % red
            
            
            n_trials_right_tmp = sum(trial_selection & R_hit_allSession==1);
            n_trials_left_tmp = sum(trial_selection & L_hit_allSession==1);
            
            subplot(5,2,n_plot); hold on
            [y x] = hist(lick_times_tmp(i_lick_right_tmp,2), 0:.1:5);
            plot(x,y/.1/n_trials_right_tmp,'b')
            [y x] = hist(lick_times_tmp(i_lick_left_tmp,2), 0:.1:5);
            plot(x,y/.1/n_trials_left_tmp,'r')
            line([0.5723 0.5723],[0 30],'color','k')
            line([1.8727  1.8727],[0 30],'color','k')
            line([3.1732 3.1732],[0 30],'color','k')
            xlim([0 4])
            ylim([0 30])
            
            
        end
    end
    
end
subplot(5,2,1); title('All mice')


%% pooled across sample and delay
figure; hold on
for i_mice = 1:n_animals
    
    
    X_type = [];
    R_early_lick = [];
    L_early_lick = [];
    n_trials = [];
    for i_aom = [0 .1 .2 .3 .4 .5]
        
        % invert the trial type based on stimulated SC
        i_select_leftSC = find(StimSide_allSession==1);
        i_select_rightSC = find(StimSide_allSession==2);
        
        Lick_ContraIpsi_allSession(i_select_leftSC,:)  = Lick_Side_allSession(i_select_leftSC);
        Lick_ContraIpsi_allSession(i_select_rightSC,:) = -(Lick_Side_allSession(i_select_rightSC)-2)+2;
        
        % select stim side
        stim_side_selection = (StimSide_allSession==i_SC_side);
        
        
        % lick early
        if i_aom ==0
            trial_selection = (AOM_data_allSession == i_aom & StimTrials_allSession >0 & Session_Index_allSession(:,1)==i_mice);
            
            X_type(end+1,1) = size(X_type,1)+1;
            R_early_lick(end+1,1) = sum(trial_selection & AOM_data_allSession == 0 & StimTrials_allSession >0 & LickEarly_allSession == 1 & Lick_ContraIpsi_allSession==3)/sum(trial_selection & AOM_data_allSession == 0 & StimTrials_allSession == 1);
            L_early_lick(end+1,1) = sum(trial_selection & AOM_data_allSession == 0 & StimTrials_allSession >0 & LickEarly_allSession == 1 & Lick_ContraIpsi_allSession==1)/sum(trial_selection & AOM_data_allSession == 0 & StimTrials_allSession == 1);
            n_trials(end+1,1) = sum(trial_selection & AOM_data_allSession == 0 & StimTrials_allSession>0);
        else
            trial_selection = (AOM_data_allSession == i_aom & (Stim_Period_allSession==1 |Stim_Period_allSession==2) & StimTrials_allSession>0 & Session_Index_allSession(:,1)==i_mice);
            if sum(trial_selection & stim_side_selection)>0
                X_type(end+1,1) = size(X_type,1)+1;
                R_early_lick(end+1,1) = sum(trial_selection & stim_side_selection & LickEarly_allSession == 1 & Lick_ContraIpsi_allSession==3)/sum(trial_selection & stim_side_selection);
                L_early_lick(end+1,1) = sum(trial_selection & stim_side_selection & LickEarly_allSession == 1 & Lick_ContraIpsi_allSession==1)/sum(trial_selection & stim_side_selection);
                n_trials(end+1,1) = sum(trial_selection & stim_side_selection);
            end
        end
        
    end
    
    if length(n_trials)>1
        plot(X_type,R_early_lick,'-b');%,'markerfacecolor','w')
        plot(X_type,L_early_lick,'-r');%,'markerfacecolor','w')
    end
    
end
   
title('both SC hemi (sample or delay stim)')
xlabel('Power')
ylabel('Fraction of early licks')
ylim([0 1])
xlim([0 6])
legend('lick contra', 'lick ipsi')

