%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 4d ALM activity in control and SC photoinhibition trials projected onto the coding dimension


clear all
close all

addpath('../func/')

load Figure_4_d_data


%% average trajectories along CD across session

% Ipsi SC stim

proj_R_training_mode1 = [];
proj_L_training_mode1 = [];
proj_R_nonStim_mode1 = [];
proj_L_nonStim_mode1 = [];
proj_R_stim1_mode1 = [];
proj_L_stim1_mode1 = [];
proj_R_stim2_mode1 = [];
proj_L_stim2_mode1 = [];

n_session = 0;
for i_session = find(RecordingSide_all(:,1)==2)'           % Ipsi SC stim
    
    n_session = n_session+1;
    
    w_CD_mode = w_CD_Ramp_all{i_session}(:,1);
    w_Ramp_mode = w_CD_Ramp_all{i_session}(:,2);
    
    % compute other modes from SVD based on actiivty variance,
    % orthogonalize w/ CD mode
    i_timebin = find(time_bins>-2.6 & time_bins<=-.4);     %only use the LDA before cue
    
    activity_matrix = activity_matrix_allSession{i_session};
    
    R = mean(activity_matrix(:,i_timebin,i_yes_trainingCorrectTrials_allSession{i_session}),3);
    L = mean(activity_matrix(:,i_timebin,i_no_trainingCorrectTrials_allSession{i_session}),3);
    activityRL = [R L];
    [u s v] = svd(activityRL');
    orthonormal_basis = Gram_Schmidt_process([w_CD_mode w_Ramp_mode v]);
    
    
    
    % project data onto the modes
    activity_matrix = activity_matrix_allSession{i_session};
    
    % training data (correct & error) (used to compute offsets)
    R = mean(activity_matrix(:,:,i_yes_training_allSession{i_session}),3);
    L = mean(activity_matrix(:,:,i_no_training_allSession{i_session}),3);
    proj_R_training = R'*orthonormal_basis;
    proj_L_training = L'*orthonormal_basis;
    
    % test non-stim data
    R = mean(activity_matrix(:,:,i_yes_nonstim_allSession{i_session}),3);
    L = mean(activity_matrix(:,:,i_no_nonstim_allSession{i_session}),3);
    proj_R_nonStim = R'*orthonormal_basis;
    proj_L_nonStim = L'*orthonormal_basis;
    
    % test sample-stim data
    R = mean(activity_matrix(:,:,i_yes_stim1_allSession{i_session}),3);
    L = mean(activity_matrix(:,:,i_no_stim1_allSession{i_session}),3);
    proj_R_stim1 = R'*orthonormal_basis;
    proj_L_stim1 = L'*orthonormal_basis;
    
    % test delay-stim data
    R = mean(activity_matrix(:,:,i_yes_stim2_allSession{i_session}),3);
    L = mean(activity_matrix(:,:,i_no_stim2_allSession{i_session}),3);
    proj_R_stim2 = R'*orthonormal_basis;
    proj_L_stim2 = L'*orthonormal_basis;
    
    
    % [1- contralateral to stimulated SC; 2- ipsi;   1-left ALM; 2-right ALM]
    
    if RecordingSide_all(i_session,2)== 1       % left ALM recording
        proj_R_training_mode1(n_session,:) = proj_R_training(:,1);
        proj_L_training_mode1(n_session,:) = proj_L_training(:,1);
        proj_R_nonStim_mode1(n_session,:) = proj_R_nonStim(:,1);
        proj_L_nonStim_mode1(n_session,:) = proj_L_nonStim(:,1);
        proj_R_stim1_mode1(n_session,:) = proj_R_stim1(:,1);
        proj_L_stim1_mode1(n_session,:) = proj_L_stim1(:,1);
        proj_R_stim2_mode1(n_session,:) = proj_R_stim2(:,1);
        proj_L_stim2_mode1(n_session,:) = proj_L_stim2(:,1);
        
        proj_R_training_mode2(n_session,:) = proj_R_training(:,2);
        proj_L_training_mode2(n_session,:) = proj_L_training(:,2);
        proj_R_nonStim_mode2(n_session,:) = proj_R_nonStim(:,2);
        proj_L_nonStim_mode2(n_session,:) = proj_L_nonStim(:,2);
        proj_R_stim1_mode2(n_session,:) = proj_R_stim1(:,2);
        proj_L_stim1_mode2(n_session,:) = proj_L_stim1(:,2);
        proj_R_stim2_mode2(n_session,:) = proj_R_stim2(:,2);
        proj_L_stim2_mode2(n_session,:) = proj_L_stim2(:,2);
        
        proj_R_training_mode3(n_session,:) = proj_R_training(:,3);
        proj_L_training_mode3(n_session,:) = proj_L_training(:,3);
        proj_R_nonStim_mode3(n_session,:) = proj_R_nonStim(:,3);
        proj_L_nonStim_mode3(n_session,:) = proj_L_nonStim(:,3);
        proj_R_stim1_mode3(n_session,:) = proj_R_stim1(:,3);
        proj_L_stim1_mode3(n_session,:) = proj_L_stim1(:,3);
        proj_R_stim2_mode3(n_session,:) = proj_R_stim2(:,3);
        proj_L_stim2_mode3(n_session,:) = proj_L_stim2(:,3);
        
    else  % righ ALM recording
        
        % flip left and right, so contra is always blue
        proj_R_training_mode1(n_session,:) = proj_L_training(:,1);
        proj_L_training_mode1(n_session,:) = proj_R_training(:,1);
        proj_R_nonStim_mode1(n_session,:) = proj_L_nonStim(:,1);
        proj_L_nonStim_mode1(n_session,:) = proj_R_nonStim(:,1);
        proj_R_stim1_mode1(n_session,:) = proj_L_stim1(:,1);
        proj_L_stim1_mode1(n_session,:) = proj_R_stim1(:,1);
        proj_R_stim2_mode1(n_session,:) = proj_L_stim2(:,1);
        proj_L_stim2_mode1(n_session,:) = proj_R_stim2(:,1);

        proj_R_training_mode2(n_session,:) = proj_R_training(:,2);
        proj_L_training_mode2(n_session,:) = proj_L_training(:,2);
        proj_R_nonStim_mode2(n_session,:) = proj_R_nonStim(:,2);
        proj_L_nonStim_mode2(n_session,:) = proj_L_nonStim(:,2);
        proj_R_stim1_mode2(n_session,:) = proj_R_stim1(:,2);
        proj_L_stim1_mode2(n_session,:) = proj_L_stim1(:,2);
        proj_R_stim2_mode2(n_session,:) = proj_R_stim2(:,2);
        proj_L_stim2_mode2(n_session,:) = proj_L_stim2(:,2);
        
        proj_R_training_mode3(n_session,:) = proj_R_training(:,3);
        proj_L_training_mode3(n_session,:) = proj_L_training(:,3);
        proj_R_nonStim_mode3(n_session,:) = proj_R_nonStim(:,3);
        proj_L_nonStim_mode3(n_session,:) = proj_L_nonStim(:,3);
        proj_R_stim1_mode3(n_session,:) = proj_R_stim1(:,3);
        proj_L_stim1_mode3(n_session,:) = proj_L_stim1(:,3);
        proj_R_stim2_mode3(n_session,:) = proj_R_stim2(:,3);
        proj_L_stim2_mode3(n_session,:) = proj_L_stim2(:,3);
        
    end


end


% CD mode
proj_training_mode1 = (proj_R_training_mode1+proj_L_training_mode1)/2;
offset_mode1 = repmat(mean(proj_training_mode1,2),1,size(proj_training_mode1,2));

proj_R_nonStim_mode1 = proj_R_nonStim_mode1-offset_mode1;
proj_L_nonStim_mode1 = proj_L_nonStim_mode1-offset_mode1;
proj_R_stim1_mode1 = proj_R_stim1_mode1-offset_mode1;
proj_L_stim1_mode1 = proj_L_stim1_mode1-offset_mode1;
proj_R_stim2_mode1 = proj_R_stim2_mode1-offset_mode1;
proj_L_stim2_mode1 = proj_L_stim2_mode1-offset_mode1;


y_max = max([mean(proj_R_nonStim_mode1) mean(proj_L_nonStim_mode1) mean(proj_R_stim1_mode1) mean(proj_L_stim1_mode1) mean(proj_R_stim2_mode1) mean(proj_L_stim2_mode1)]);
y_min = min([mean(proj_R_nonStim_mode1) mean(proj_L_nonStim_mode1) mean(proj_R_stim1_mode1) mean(proj_L_stim1_mode1) mean(proj_R_stim2_mode1) mean(proj_L_stim2_mode1)]);
y_max = y_max+abs(y_max)*.5;
y_min = y_min-abs(y_min)*.5;


figure
subplot(1,2,1);
func_plot_mean_and_sem(time_bins+.2,proj_R_nonStim_mode1, 'b', [.7 .7 1], 'n')
func_plot_mean_and_sem(time_bins+.2,proj_L_nonStim_mode1, 'r', [1 .7 .7], 'n')
line([-2.6 -2.6],[y_min y_max],'color','k')
line([-1.3 -1.3],[y_min y_max],'color','k')
line([0 0],[y_min y_max],'color','k')
xlim([-2.8 1])
ylim([y_min y_max])
title('Ipsi SC photoinhibition, CD mode')
legend(' ','contraLick',' ','ipsiLick')


subplot(1,2,2);
func_plot_mean_and_sem(time_bins+.2,proj_R_stim2_mode1, 'b', [.7 .7 1], 'n')
func_plot_mean_and_sem(time_bins+.2,proj_L_stim2_mode1, 'r', [1 .7 .7], 'n')
plot(time_bins+.2,mean(proj_R_nonStim_mode1), 'b')
plot(time_bins+.2,mean(proj_L_nonStim_mode1), 'r')
line([-2.6 -2.6],[y_min y_max],'color','k')
line([-1.3 -1.3],[y_min y_max],'color','k')
line([0 0],[y_min y_max],'color','k')
line([-1.3 0],[y_max y_max]*.8,'color','c','linewidth',3)
xlim([-2.8 1])
ylim([y_min y_max])




% Ramp mode
proj_training_mode2 = (proj_R_training_mode2+proj_L_training_mode2)/2;
offset_mode2 = repmat(mean(proj_training_mode2,2),1,size(proj_training_mode2,2));

proj_R_nonStim_mode2 = proj_R_nonStim_mode2-offset_mode2;
proj_L_nonStim_mode2 = proj_L_nonStim_mode2-offset_mode2;
proj_R_stim1_mode2 = proj_R_stim1_mode2-offset_mode2;
proj_L_stim1_mode2 = proj_L_stim1_mode2-offset_mode2;
proj_R_stim2_mode2 = proj_R_stim2_mode2-offset_mode2;
proj_L_stim2_mode2 = proj_L_stim2_mode2-offset_mode2;


y_max = max([mean(proj_R_nonStim_mode2) mean(proj_L_nonStim_mode2) mean(proj_R_stim1_mode2) mean(proj_L_stim1_mode2) mean(proj_R_stim2_mode2) mean(proj_L_stim2_mode2)]);
y_min = min([mean(proj_R_nonStim_mode2) mean(proj_L_nonStim_mode2) mean(proj_R_stim1_mode2) mean(proj_L_stim1_mode2) mean(proj_R_stim2_mode2) mean(proj_L_stim2_mode2)]);
y_max = y_max+abs(y_max)*.5;
y_min = y_min-abs(y_min)*.5;


figure
subplot(1,2,1);
func_plot_mean_and_sem(time_bins+.2,proj_R_nonStim_mode2, 'b', [.7 .7 1], 'n')
func_plot_mean_and_sem(time_bins+.2,proj_L_nonStim_mode2, 'r', [1 .7 .7], 'n')
line([-2.6 -2.6],[y_min y_max],'color','k')
line([-1.3 -1.3],[y_min y_max],'color','k')
line([0 0],[y_min y_max],'color','k')
xlim([-2.8 1])
ylim([y_min y_max])
title('Ipsi SC photoinhibition, Ramp mode')
legend(' ','contraLick',' ','ipsiLick')

subplot(1,2,2);
func_plot_mean_and_sem(time_bins+.2,proj_R_stim2_mode2, 'b', [.7 .7 1], 'n')
func_plot_mean_and_sem(time_bins+.2,proj_L_stim2_mode2, 'r', [1 .7 .7], 'n')
plot(time_bins+.2,mean(proj_R_nonStim_mode2), 'b')
plot(time_bins+.2,mean(proj_L_nonStim_mode2), 'r')
line([-2.6 -2.6],[y_min y_max],'color','k')
line([-1.3 -1.3],[y_min y_max],'color','k')
line([0 0],[y_min y_max],'color','k')
line([-1.3 0],[y_max y_max]*.8,'color','c','linewidth',3)
xlim([-2.8 1])
ylim([y_min y_max])



% First PC mode
proj_training_mode3 = (proj_R_training_mode3+proj_L_training_mode3)/2;
offset_mode3 = repmat(mean(proj_training_mode3,2),1,size(proj_training_mode3,2));

proj_R_nonStim_mode3 = proj_R_nonStim_mode3-offset_mode3;
proj_L_nonStim_mode3 = proj_L_nonStim_mode3-offset_mode3;
proj_R_stim1_mode3 = proj_R_stim1_mode3-offset_mode3;
proj_L_stim1_mode3 = proj_L_stim1_mode3-offset_mode3;
proj_R_stim2_mode3 = proj_R_stim2_mode3-offset_mode3;
proj_L_stim2_mode3 = proj_L_stim2_mode3-offset_mode3;


y_max = max([mean(proj_R_nonStim_mode3) mean(proj_L_nonStim_mode3) mean(proj_R_stim1_mode3) mean(proj_L_stim1_mode3) mean(proj_R_stim2_mode3) mean(proj_L_stim2_mode3)]);
y_min = min([mean(proj_R_nonStim_mode3) mean(proj_L_nonStim_mode3) mean(proj_R_stim1_mode3) mean(proj_L_stim1_mode3) mean(proj_R_stim2_mode3) mean(proj_L_stim2_mode3)]);
y_max = y_max+abs(y_max)*.5;
y_min = y_min-abs(y_min)*.5;


figure
subplot(1,2,1);
func_plot_mean_and_sem(time_bins+.2,proj_R_nonStim_mode3, 'b', [.7 .7 1], 'n')
func_plot_mean_and_sem(time_bins+.2,proj_L_nonStim_mode3, 'r', [1 .7 .7], 'n')
line([-2.6 -2.6],[y_min y_max],'color','k')
line([-1.3 -1.3],[y_min y_max],'color','k')
line([0 0],[y_min y_max],'color','k')
xlim([-2.8 1])
ylim([y_min y_max])
title('Ipsi SC photoinhibition, 1st PC mode')
legend(' ','contraLick',' ','ipsiLick')


subplot(1,2,2);
func_plot_mean_and_sem(time_bins+.2,proj_R_stim2_mode3, 'b', [.7 .7 1], 'n')
func_plot_mean_and_sem(time_bins+.2,proj_L_stim2_mode3, 'r', [1 .7 .7], 'n')
plot(time_bins+.2,mean(proj_R_nonStim_mode3), 'b')
plot(time_bins+.2,mean(proj_L_nonStim_mode3), 'r')
line([-2.6 -2.6],[y_min y_max],'color','k')
line([-1.3 -1.3],[y_min y_max],'color','k')
line([0 0],[y_min y_max],'color','k')
line([-1.3 0],[y_max y_max]*.8,'color','c','linewidth',3)
xlim([-2.8 1])
ylim([y_min y_max])


