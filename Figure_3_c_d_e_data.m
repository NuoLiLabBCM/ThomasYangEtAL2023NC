clear all
close all

load Figure_3_c_d_e_data



%% Group neurons into Lick Right preferring and Lick Left preferring, examine their single trial dynamics together

corr_r_yes = [];
corr_r_no = [];
R_neurons_singleTrials_yes_allSesssions = [];
L_neurons_singleTrials_yes_allSesssions = [];
R_neurons_singleTrials_no_allSesssions = [];
L_neurons_singleTrials_no_allSesssions = [];
session_ID_yes = [];
session_ID_no = [];

R_neurons_yes_allSesssions = [];
L_neurons_yes_allSesssions = [];
R_neurons_no_allSesssions = [];
L_neurons_no_allSesssions = [];

i_session_selected = [];

figure
n_R_list = [];
n_L_list = [];
for i_session = 1:size(activity_matrix_allSession,1)
        
    i_timebin = find(time_bins>-3.2 & time_bins<=-.4);     %only use the LDA before cue
    
    activity_matrix = activity_matrix_allSession{i_session};
    
    
    spk_count_yes_screen = spk_count_yes_screen_allSession{i_session,1};    % neuron x 5, [whole_trial   sample   delay   resposne]
    spk_count_no_screen = spk_count_no_screen_allSession{i_session,1};
    
    sig_selective = sig_selective_allSession{i_session,1};

    i_yes_screen_trial = i_yes_screen_allSession{i_session,1};
    i_no_screen_trial = i_no_screen_allSession{i_session,1};
    i_yes_correct_trial = i_yes_correct_allSession{i_session,1};  % indices for non-stim trials, independent test data only
    i_no_correct_trial = i_no_correct_allSession{i_session,1};
    i_yes_error_trial = i_yes_error_allSession{i_session,1};
    i_no_error_trial = i_no_error_allSession{i_session,1};

    
    %average_activity = nanmean(activity_matrix,3);      % average spike rate across all trials types, this will be substracted at each time point (detrending the data) to examine competition dynamics more closely
    %average_activity = nanmean(activity_matrix(:,:,[i_yes_screen_trial; i_no_screen_trial; i_yes_error_trial; i_no_error_trial]),3);      % average spike rate across all trials types, this will be substracted at each time point (detrending the data) to examine competition dynamics more closely
    average_activity = nanmean(activity_matrix(:,:,[i_yes_screen_trial i_no_screen_trial]),3);      % average spike rate across all trials types, this will be substracted at each time point (detrending the data) to examine competition dynamics more closely
    
    sig_selective = sig_selective<.01;      % selective neurons only

    i_selR_unit = find(spk_count_yes_screen(:,3)>spk_count_no_screen(:,3) & sum(sig_selective(:,1:3),2)>0);     % yes perferring neurons,  this is not yet flipped to contra/ipsi
    i_selL_unit = find(spk_count_yes_screen(:,3)<spk_count_no_screen(:,3) & sum(sig_selective(:,1:3),2)>0);     % no perferring neurons

    if size(i_selR_unit,1)>=5 & size(i_selL_unit,1)>=5
        
        R_neurons_yes_iSesssions = [];
        L_neurons_yes_iSesssions = [];
        R_neurons_no_iSesssions = [];
        L_neurons_no_iSesssions = [];
        
        n_plot = 0;
        for i_trial = i_yes_correct_trial'
            n_plot = n_plot+1;
            
            % mean activity of yes and no perferring neurons
            R_neurons = mean(activity_matrix(i_selR_unit,i_timebin,i_trial)-average_activity(i_selR_unit,i_timebin));
            L_neurons = mean(activity_matrix(i_selL_unit,i_timebin,i_trial)-average_activity(i_selL_unit,i_timebin));
            
            corr_r_yes(end+1,1) = corr(R_neurons',L_neurons');
            
            R_neurons_yes_iSesssions(end+1,:) = R_neurons;
            L_neurons_yes_iSesssions(end+1,:) = L_neurons;
            
        end
        
        for i_trial = i_no_correct_trial'
            n_plot = n_plot+1;
            
            % mean activity of yes and no perferring neurons
            R_neurons = mean(activity_matrix(i_selR_unit,i_timebin,i_trial)-average_activity(i_selR_unit,i_timebin));
            L_neurons = mean(activity_matrix(i_selL_unit,i_timebin,i_trial)-average_activity(i_selL_unit,i_timebin));
            
            corr_r_no(end+1,1) = corr(R_neurons',L_neurons');

            R_neurons_no_iSesssions(end+1,:) = R_neurons;
            L_neurons_no_iSesssions(end+1,:) = L_neurons;
        end
        
        R_neurons_singleTrials_yes_allSesssions = cat(1,R_neurons_singleTrials_yes_allSesssions, R_neurons_yes_iSesssions);
        L_neurons_singleTrials_yes_allSesssions = cat(1,L_neurons_singleTrials_yes_allSesssions, L_neurons_yes_iSesssions);
        R_neurons_singleTrials_no_allSesssions = cat(1,R_neurons_singleTrials_no_allSesssions, R_neurons_no_iSesssions);
        L_neurons_singleTrials_no_allSesssions = cat(1,L_neurons_singleTrials_no_allSesssions, L_neurons_no_iSesssions);
        session_ID_yes = cat(1,session_ID_yes,ones(size(R_neurons_yes_iSesssions,1),1)*i_session);
        session_ID_no = cat(1,session_ID_no,ones(size(L_neurons_no_iSesssions,1),1)*i_session);
        
        
        population_sel1 = abs(mean(R_neurons_yes_iSesssions(:,end))-mean(R_neurons_no_iSesssions(:,end)));
        population_sel2 = abs(mean(L_neurons_yes_iSesssions(:,end))-mean(L_neurons_no_iSesssions(:,end)));
        
        if population_sel1>2 & population_sel2>2
            i_session_selected = [i_session_selected; i_session];
            n_R_list(end+1,:) = size(i_selR_unit,1);
            n_L_list(end+1,:) = size(i_selL_unit,1);
            
            hold on
            plot(mean(R_neurons_yes_iSesssions),mean(L_neurons_yes_iSesssions),'b');
            plot(mean(R_neurons_yes_iSesssions(:,end)),mean(L_neurons_yes_iSesssions(:,end)),'ow','markerfacecolor','b');
            plot(mean(R_neurons_no_iSesssions),mean(L_neurons_no_iSesssions),'r');
            plot(mean(R_neurons_no_iSesssions(:,end)),mean(L_neurons_no_iSesssions(:,end)),'ow','markerfacecolor','r');
        end
    end
    
end
title('All sessions')
xlabel('Lick right pref. neurons')
ylabel('Lick left pref. neurons')





% i_session = 8;
% i_session = 9;
i_session = 11;%21;

% lick right trials
figure;
R_neurons = R_neurons_singleTrials_yes_allSesssions(session_ID_yes==i_session,:);
L_neurons = L_neurons_singleTrials_yes_allSesssions(session_ID_yes==i_session,:);
corr_r_yes_tmp = corr_r_yes(session_ID_yes==i_session,:);
[dummy i_sort] = sort(corr_r_yes_tmp);
n_plot = 0;
for i_trial = i_sort(1:5)'
    n_plot = n_plot+1;
    subplot(4,5,n_plot); hold on
    plot(time_bins(i_timebin),R_neurons(i_trial,:),'b');
    plot(time_bins(i_timebin),L_neurons(i_trial,:),'r');

    subplot(4,5,n_plot+5); hold on
    plot(R_neurons(i_trial,:), L_neurons(i_trial,:),'b');
    plot(R_neurons(i_trial,end), L_neurons(i_trial,end),'ow','markerfacecolor','b');
end

% lick left trials
R_neurons = R_neurons_singleTrials_no_allSesssions(session_ID_no==i_session,:);
L_neurons = L_neurons_singleTrials_no_allSesssions(session_ID_no==i_session,:);
corr_r_no_tmp = corr_r_no(session_ID_no==i_session,:);
[dummy i_sort] = sort(corr_r_no_tmp);
n_plot = 0;
for i_trial = i_sort(1:5)'
    n_plot = n_plot+1;
    subplot(4,5,n_plot+10); hold on
    plot(time_bins(i_timebin),R_neurons(i_trial,:),'b');
    plot(time_bins(i_timebin),L_neurons(i_trial,:),'r');

    subplot(4,5,n_plot+15); hold on
    plot(R_neurons(i_trial,:), L_neurons(i_trial,:),'r');
    plot(R_neurons(i_trial,end), L_neurons(i_trial,end),'ow','markerfacecolor','r');
end

subplot(4,5,1);
title('Example Lick Right trials')
xlabel('Time (s)')
ylabel('delta FR (spk/s)')

subplot(4,5,5)
legend('Lick right pref neurons','Lick left pref neurons')

subplot(4,5,6);
xlabel('Right pref neurons')
ylabel('Left pref neurons')

subplot(4,5,11);
title('Example Lick Left trials')




figure(10); hold on
[y x] = hist(corr_r_yes(ismember(session_ID_yes,i_session_selected)));
plot(x,y/sum(y),'b');
[y x] = hist(corr_r_no(ismember(session_ID_no,i_session_selected)));
plot(x,y/sum(y),'r');
xlabel('Correlation between Left and Right perf neurons')
ylabel('Number of trials')

disp(['======================']);
disp(['yes: ',num2str(mean(corr_r_yes(ismember(session_ID_yes,i_session_selected)))),'+',num2str(std(corr_r_yes(ismember(session_ID_yes,i_session_selected)))/sqrt(size(corr_r_yes(ismember(session_ID_yes,i_session_selected)),1)))]);
disp(['no: ',num2str(mean(corr_r_no(ismember(session_ID_no,i_session_selected)))),'+',num2str(std(corr_r_no(ismember(session_ID_no,i_session_selected)))/sqrt(size(corr_r_no(ismember(session_ID_no,i_session_selected)),1)))]);



%% Shuffle control, random groups of neurons

corr_r_yes = [];
corr_r_no = [];
R_neurons_singleTrials_yes_allSesssions = [];
L_neurons_singleTrials_yes_allSesssions = [];
R_neurons_singleTrials_no_allSesssions = [];
L_neurons_singleTrials_no_allSesssions = [];
session_ID_yes = [];
session_ID_no = [];

R_neurons_yes_allSesssions = [];
L_neurons_yes_allSesssions = [];
R_neurons_no_allSesssions = [];
L_neurons_no_allSesssions = [];

figure
for i_session = i_session_selected'%1:size(activity_matrix_allSession,1)
        
    i_timebin = find(time_bins>-3.2 & time_bins<=-.4);     %only use the LDA before cue
    
    activity_matrix = activity_matrix_allSession{i_session};
    
    
    spk_count_yes_screen = spk_count_yes_screen_allSession{i_session,1};    % neuron x 5, [whole_trial   sample   delay   resposne]
    spk_count_no_screen = spk_count_no_screen_allSession{i_session,1};
    
    sig_selective = sig_selective_allSession{i_session,1};

    
    i_yes_screen_trial = i_yes_screen_allSession{i_session,1};
    i_no_screen_trial = i_no_screen_allSession{i_session,1};
    i_yes_correct_trial = i_yes_correct_allSession{i_session,1};  % indices for non-stim trials, independent test data only
    i_no_correct_trial = i_no_correct_allSession{i_session,1};
    i_yes_error_trial = i_yes_error_allSession{i_session,1};
    i_no_error_trial = i_no_error_allSession{i_session,1};

    
    %average_activity = nanmean(activity_matrix,3);      % average spike rate across all trials types, this will be substracted at each time point (detrending the data) to examine competition dynamics more closely
    %average_activity = nanmean(activity_matrix(:,:,[i_yes_screen_trial; i_no_screen_trial; i_yes_error_trial; i_no_error_trial]),3);      % average spike rate across all trials types, this will be substracted at each time point (detrending the data) to examine competition dynamics more closely
    average_activity = nanmean(activity_matrix(:,:,[i_yes_screen_trial i_no_screen_trial]),3);      % average spike rate across all trials types, this will be substracted at each time point (detrending the data) to examine competition dynamics more closely
    
    sig_selective = sig_selective<.01;      % selective neurons only

    %i_selR_unit = find(spk_count_yes_screen(:,3)>spk_count_no_screen(:,3) & sum(sig_selective(:,1:3),2)>0);     % yes perferring neurons,  this is not yet flipped to contra/ipsi
    %i_selL_unit = find(spk_count_yes_screen(:,3)<spk_count_no_screen(:,3) & sum(sig_selective(:,1:3),2)>0);     % no perferring neurons
    i_sell_unit = find(sum(sig_selective(:,1:3),2)>0); 
    i_selR_unit = i_sell_unit(1:round(length(i_sell_unit)/2));     % yes perferring neurons,  this is not yet flipped to contra/ipsi
    i_selL_unit = i_sell_unit(round(length(i_sell_unit)/2)+1:end);     % no perferring neurons

    if size(i_selR_unit,1)>=5 & size(i_selL_unit,1)>=5
        
        R_neurons_yes_iSesssions = [];
        L_neurons_yes_iSesssions = [];
        R_neurons_no_iSesssions = [];
        L_neurons_no_iSesssions = [];
        
        n_plot = 0;
        for i_trial = i_yes_correct_trial'
            n_plot = n_plot+1;
            
            % mean activity of yes and no perferring neurons
            R_neurons = mean(activity_matrix(i_selR_unit,i_timebin,i_trial)-average_activity(i_selR_unit,i_timebin));
            L_neurons = mean(activity_matrix(i_selL_unit,i_timebin,i_trial)-average_activity(i_selL_unit,i_timebin));
            
            corr_r_yes(end+1,1) = corr(R_neurons',L_neurons');
            
            R_neurons_yes_iSesssions(end+1,:) = R_neurons;
            L_neurons_yes_iSesssions(end+1,:) = L_neurons;
            
        end
        
        for i_trial = i_no_correct_trial'
            n_plot = n_plot+1;
            
            % mean activity of yes and no perferring neurons
            R_neurons = mean(activity_matrix(i_selR_unit,i_timebin,i_trial)-average_activity(i_selR_unit,i_timebin));
            L_neurons = mean(activity_matrix(i_selL_unit,i_timebin,i_trial)-average_activity(i_selL_unit,i_timebin));
            
            corr_r_no(end+1,1) = corr(R_neurons',L_neurons');

            R_neurons_no_iSesssions(end+1,:) = R_neurons;
            L_neurons_no_iSesssions(end+1,:) = L_neurons;
        end
        
        R_neurons_singleTrials_yes_allSesssions = cat(1,R_neurons_singleTrials_yes_allSesssions, R_neurons_yes_iSesssions);
        L_neurons_singleTrials_yes_allSesssions = cat(1,L_neurons_singleTrials_yes_allSesssions, L_neurons_yes_iSesssions);
        R_neurons_singleTrials_no_allSesssions = cat(1,R_neurons_singleTrials_no_allSesssions, R_neurons_no_iSesssions);
        L_neurons_singleTrials_no_allSesssions = cat(1,L_neurons_singleTrials_no_allSesssions, L_neurons_no_iSesssions);
        session_ID_yes = cat(1,session_ID_yes,ones(size(R_neurons_yes_iSesssions,1),1)*i_session);
        session_ID_no = cat(1,session_ID_no,ones(size(L_neurons_no_iSesssions,1),1)*i_session);
        
    end
    
    hold on
    plot(mean(R_neurons_yes_iSesssions),mean(L_neurons_yes_iSesssions),'b');
    plot(mean(R_neurons_yes_iSesssions(:,end)),mean(L_neurons_yes_iSesssions(:,end)),'ow','markerfacecolor','b');
    plot(mean(R_neurons_no_iSesssions),mean(L_neurons_no_iSesssions),'r');
    plot(mean(R_neurons_no_iSesssions(:,end)),mean(L_neurons_no_iSesssions(:,end)),'ow','markerfacecolor','r');
    
end
title('All sessions, random neuron group')
xlabel('Neurons group 1')
ylabel('Neurons group 2')


figure(10); hold on
[y x] = hist([corr_r_yes; corr_r_no]);
plot(x,y/sum(y),'k');
xlabel('Correlation between Left and Right perf neurons')
ylabel('Number of trials')
legend('Lick right trials, R vs. L pop.','Lick right trials, R vs. L pop.','Null, shuffl neuron grouping')


disp(['shuffle: ',num2str(mean([corr_r_yes; corr_r_no])),'+',num2str(std([corr_r_yes; corr_r_no])/sqrt(size([corr_r_yes; corr_r_no],1)))]);

