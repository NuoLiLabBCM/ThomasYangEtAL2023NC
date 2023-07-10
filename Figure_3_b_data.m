clear all
close all

load Figure-3_b_data


sig_selective_p(isnan(sig_selective_p))=1;

sig_selective = sig_selective_p<0.01;%sig_selective_p<(0.05/size(sig_selective_p,1)/3); % <----------- This p value can be adjusted



%% ---------------- flip all units into the left hemisphere ------------
i_units = find(RecordingSide_all==2);

spk_count_no_tmp = spk_count_yes_screen(i_units,:);
spk_count_yes_tmp = spk_count_no_screen(i_units,:);
spk_count_yes_screen(i_units,:) = spk_count_yes_tmp;
spk_count_no_screen(i_units,:) = spk_count_no_tmp;

spk_count_no_tmp = spk_count_yes_correct_all(i_units,:);
spk_count_yes_tmp = spk_count_no_correct_all(i_units,:);
spk_count_yes_correct_all(i_units,:) = spk_count_yes_tmp;
spk_count_no_correct_all(i_units,:) = spk_count_no_tmp;

spk_count_no_tmp = spk_count_yes_error_all(i_units,:);
spk_count_yes_tmp = spk_count_no_error_all(i_units,:);
spk_count_yes_error_all(i_units,:) = spk_count_yes_tmp;
spk_count_no_error_all(i_units,:) = spk_count_no_tmp;

PSTH_no_tmp = PSTH_yes_test_cue_aligned(i_units,:);
PSTH_yes_tmp = PSTH_no_test_cue_aligned(i_units,:);
PSTH_yes_test_cue_aligned(i_units,:) = PSTH_yes_tmp;
PSTH_no_test_cue_aligned(i_units,:) = PSTH_no_tmp;


PSTH_no_tmp = PSTH_yes_correct_cue_aligned(i_units,:);
PSTH_yes_tmp = PSTH_no_correct_cue_aligned(i_units,:);
PSTH_yes_correct_cue_aligned(i_units,:) = PSTH_yes_tmp;
PSTH_no_correct_cue_aligned(i_units,:) = PSTH_no_tmp;

PSTH_no_tmp = PSTH_yes_error_cue_aligned(i_units,:);
PSTH_yes_tmp = PSTH_no_error_cue_aligned(i_units,:);
PSTH_yes_error_cue_aligned(i_units,:) = PSTH_yes_tmp;
PSTH_no_error_cue_aligned(i_units,:) = PSTH_no_tmp;


spike_times_tmp1 = spk_times_yes_correct_all(i_units,1);
spike_times_tmp2 = spk_times_no_correct_all(i_units,1);
spike_times_tmp3 = spk_times_yes_error_all(i_units,1);
spike_times_tmp4 = spk_times_no_error_all(i_units,1);

spk_times_no_correct_all(i_units,1) = spike_times_tmp1;
spk_times_yes_correct_all(i_units,1) = spike_times_tmp2;
spk_times_no_error_all(i_units,1) = spike_times_tmp3;
spk_times_yes_error_all(i_units,1) = spike_times_tmp4;




%% selectivity, fraction of neurons selective, all of SC neurons, separately sort selectivity based on sample, delay, and response epoch
FR_pref = spk_count_yes_screen-spk_count_no_screen;         % [whole_trial  sample  delay   response]

frac_prep_only = sum((sig_selective(:,1)|sig_selective(:,2)) & ~(sig_selective(:,3)))/size(sig_selective,1);
frac_prep_mov = sum((sig_selective(:,1)|sig_selective(:,2)) & (sig_selective(:,3)))/size(sig_selective,1);
frac_mov_only = sum(~(sig_selective(:,1)|sig_selective(:,2)) & (sig_selective(:,3)))/size(sig_selective,1);

figure;
subplot(1,4,1); hold on
bar(1,frac_prep_only,'w');
bar(2,frac_prep_mov,'faceColor',[.6 .6 .6]);
bar(3,frac_mov_only,'k');
legend('Type 1','Type 2','Type 3');
title('All neurons in ALM')
ylim([0 .25])

for i_epoch = 1:3       % selectivity, preferred-nonpreferred, determined based on each epoch separately
    
    i_selective = sig_selective(:,i_epoch); % find cells that are selectivity during specific epoch
    
    i_n_trial =  min(N_trials_all(:,[1 2]),[],2)>=20;                        % find cells that have more than 10 trials in each trial type
    
    PSTH_pref = [];
    PSTH_nonpref = [];
    for i_cell = find(i_selective & i_n_trial)'
        
        if FR_pref(i_cell,1+i_epoch)>0
            
            PSTH_pref(end+1,:)            = PSTH_yes_test_cue_aligned(i_cell,:);
            PSTH_nonpref(end+1,:)         = PSTH_no_test_cue_aligned(i_cell,:);
            
            
        elseif FR_pref(i_cell,1+i_epoch)<0
            
            PSTH_pref(end+1,:)              = PSTH_no_test_cue_aligned(i_cell,:);
            PSTH_nonpref(end+1,:)           = PSTH_yes_test_cue_aligned(i_cell,:);
            
        end
    end
    
    selectivity = PSTH_pref-PSTH_nonpref;
    
    subplot(1,4,1+i_epoch); hold on
    func_plot_mean_and_sem(t,selectivity, 'k', [.7 .7 .7], 'n');
    line([-3.5 2],[0 0],'color','k')
    line([0 0],[-1 8],'color','k')
    line([-1.3 -1.3],[-1 8],'color','k')
    line([-2.6 -2.6],[-1 8],'color','k')
    xlabel('time (s)')
    ylabel('selectivity (spk/s)')
    xlim([-3.5 2])
    ylim([-1 8])
    if i_epoch == 1
        title('Based on Sample selectivity')
    elseif i_epoch == 2
        title('Delay selectivity')
    elseif i_epoch == 3
        title('Response selectivity')
    end
    
end



%% selectivity, fraction of contra vs. ipsi, all of SC neurons, separately sort selectivity based on sample, delay, and response epoch
FR_pref = spk_count_yes_screen-spk_count_no_screen;

figure;
for i_epoch = 1:3    % separately sort selectivity based on sample, delay, and response epoch
    
    i_selective = sig_selective(:,i_epoch); % find cells that are selectivity during that epoch
    
    frac_contra = sum(i_selective & FR_pref(:,1+i_epoch)>0)/sum(i_selective);
    frac_ipsi = sum(i_selective & FR_pref(:,1+i_epoch)<0)/sum(i_selective);

    
    frac_contra_iBtstrp = [];
    frac_ipsi_iBtstrp = [];
    for i_btstrp = 1:10000
        
        % resample
        Mice_ID = unique(Mice_all(:,1));
        Mice_ID_iBtstrp = randsample(Mice_ID,length(Mice_ID),'true');
        
        % build resampled dataset
        sig_selective_iBtstrp = [];
        FR_pref_iBtstrp = [];
        for i_Mice = Mice_ID_iBtstrp'
            
            sig_selective_iBtstrp = cat(1, sig_selective_iBtstrp, sig_selective(Mice_all(:,1)==i_Mice,:));
            FR_pref_iBtstrp = cat(1, FR_pref_iBtstrp, FR_pref(Mice_all(:,1)==i_Mice,:));
            
        end
        
        i_selective_btstrp = (sig_selective_iBtstrp(:,1)|sig_selective_iBtstrp(:,2)|sig_selective_iBtstrp(:,3));
        frac_contra_iBtstrp(i_btstrp,1) = sum(i_selective_btstrp & FR_pref_iBtstrp(:,1+i_epoch)>0)/sum(i_selective_btstrp);
        frac_ipsi_iBtstrp(i_btstrp,1) = sum(i_selective_btstrp & FR_pref_iBtstrp(:,1+i_epoch)<0)/sum(i_selective_btstrp);
        
    end        
    
    
    subplot(3,5,1+(i_epoch-1)*5); hold on
    bar(1,frac_contra,'b');
    bar(2,frac_ipsi,'r');
    errorbar(1, frac_contra, std(frac_contra_iBtstrp),'b');
    errorbar(2, frac_ipsi, std(frac_ipsi_iBtstrp),'r');
    ylim([0 1])
    p_value = sum(frac_contra_iBtstrp>frac_ipsi_iBtstrp)/size(frac_contra_iBtstrp,1);
    if p_value>.5
        p_value = 1-p_value;
    end
    if i_epoch == 1
        title(['All neurons in ALM, Sample, p=',num2str(p_value)])
    elseif i_epoch == 2
        title(['Delay, p=',num2str(p_value)])
    elseif i_epoch == 3
        title(['Response, p=',num2str(p_value)])
    end
    
    
    % contra & ipsi selective population
    i_contra = (i_selective & FR_pref(:,1+i_epoch)>0);
    i_ipsi = (i_selective & FR_pref(:,1+i_epoch)<0);
    
    subplot(3,5,[2 3]+(i_epoch-1)*5)
    func_plot_mean_and_sem(t,PSTH_yes_correct_cue_aligned(i_contra,:), 'b', [.7 .7 1], 'n');
    func_plot_mean_and_sem(t,PSTH_no_correct_cue_aligned(i_contra,:), 'r', [1 .7 .7], 'n');
    line([0 0],[-5 40],'color','k')
    line([-1.3 -1.3],[0 40],'color','k')
    line([-2.6 -2.6],[0 40],'color','k')
    xlabel('time (s)')
    ylabel('Spike rate (spk/s)')
    xlim([-3 1.6])
    ylim([0 20])
    title('contra preferring')
    
    subplot(3,5,[4 5]+(i_epoch-1)*5)
    func_plot_mean_and_sem(t,PSTH_yes_correct_cue_aligned(i_ipsi,:), 'b', [.7 .7 1], 'n');
    func_plot_mean_and_sem(t,PSTH_no_correct_cue_aligned(i_ipsi,:), 'r', [1 .7 .7], 'n');
    line([0 0],[0 40],'color','k')
    line([-1.3 -1.3],[-5 40],'color','k')
    line([-2.6 -2.6],[0 40],'color','k')
    xlabel('time (s)')
    ylabel('Spike rate (spk/s)')
    xlim([-3 1.6])
    ylim([0 20])
    title('ipsi preferring')
    
end




%% contra vs ipsi, whole trial or epoch
FR_pref = spk_count_yes_correct_all - spk_count_no_correct_all;


figure;

% whole trial
i_selective = (sig_selective(:,1)|sig_selective(:,2)|sig_selective(:,3));
frac_contra_SDR = sum((i_selective & FR_pref(:,1)>0 & RecordingSide_all==1)|(i_selective & FR_pref(:,1)<0 & RecordingSide_all==2))/sum((i_selective & RecordingSide_all==1)|(i_selective & RecordingSide_all==2));
frac_ipsi_SDR = sum((i_selective & FR_pref(:,1)<0 & RecordingSide_all==1)|(i_selective & FR_pref(:,1)>0 & RecordingSide_all==2))/sum((i_selective & RecordingSide_all==1)|(i_selective & RecordingSide_all==2));

% sample epoch
i_selective = (sig_selective(:,1));
frac_contra_S = sum((i_selective & FR_pref(:,2)>0 & RecordingSide_all==1)|(i_selective & FR_pref(:,2)<0 & RecordingSide_all==2))/sum((i_selective & RecordingSide_all==1)|(i_selective & RecordingSide_all==2));
frac_ipsi_S = sum((i_selective & FR_pref(:,2)<0 & RecordingSide_all==1)|(i_selective & FR_pref(:,2)>0 & RecordingSide_all==2))/sum((i_selective & RecordingSide_all==1)|(i_selective & RecordingSide_all==2));

% delay epoch
i_selective = (sig_selective(:,2));
frac_contra_D = sum((i_selective & FR_pref(:,3)>0 & RecordingSide_all==1)|(i_selective & FR_pref(:,3)<0 & RecordingSide_all==2))/sum((i_selective & RecordingSide_all==1)|(i_selective & RecordingSide_all==2));
frac_ipsi_D = sum((i_selective & FR_pref(:,3)<0 & RecordingSide_all==1)|(i_selective & FR_pref(:,3)>0 & RecordingSide_all==2))/sum((i_selective & RecordingSide_all==1)|(i_selective & RecordingSide_all==2));

% response epoch
i_selective = (sig_selective(:,3));
frac_contra_R = sum((i_selective & FR_pref(:,4)>0 & RecordingSide_all==1)|(i_selective & FR_pref(:,4)<0 & RecordingSide_all==2))/sum((i_selective & RecordingSide_all==1)|(i_selective & RecordingSide_all==2));
frac_ipsi_R = sum((i_selective & FR_pref(:,4)<0 & RecordingSide_all==1)|(i_selective & FR_pref(:,4)>0 & RecordingSide_all==2))/sum((i_selective & RecordingSide_all==1)|(i_selective & RecordingSide_all==2));



frac_contra_SDR_iBtstrp = [];
frac_contra_S_iBtstrp = [];
frac_contra_D_iBtstrp = [];
frac_contra_R_iBtstrp = [];
frac_ipsi_SDR_iBtstrp = [];
frac_ipsi_S_iBtstrp = [];
frac_ipsi_D_iBtstrp = [];
frac_ipsi_R_iBtstrp = [];
for i_btstrp = 1:10000
    
    % resample
    Mice_ID = unique(Mice_all(:,1));
    Mice_ID_iBtstrp = randsample(Mice_ID,length(Mice_ID),'true');
    
    % build resampled dataset
    sig_selective_iBtstrp = [];
    FR_pref_iBtstrp = [];
    RecordingSide_iBtstrp = [];
    for i_Mice = Mice_ID_iBtstrp'
        
        sig_selective_iBtstrp = cat(1, sig_selective_iBtstrp, sig_selective(Mice_all(:,1)==i_Mice,:));
        FR_pref_iBtstrp = cat(1, FR_pref_iBtstrp, FR_pref(Mice_all(:,1)==i_Mice,:));
        RecordingSide_iBtstrp = cat(1, RecordingSide_iBtstrp, RecordingSide_all(Mice_all(:,1)==i_Mice,:));
        
    end
    
    % whole trial
    i_selective = (sig_selective_iBtstrp(:,1)|sig_selective_iBtstrp(:,2)|sig_selective_iBtstrp(:,3));
    frac_contra_SDR_iBtstrp(i_btstrp,:) = sum((i_selective & FR_pref_iBtstrp(:,1)>0 & RecordingSide_iBtstrp==1)|(i_selective & FR_pref_iBtstrp(:,1)<0 & RecordingSide_iBtstrp==2))/sum((i_selective & RecordingSide_iBtstrp==1)|(i_selective & RecordingSide_iBtstrp==2));
    frac_ipsi_SDR_iBtstrp(i_btstrp,:) = sum((i_selective & FR_pref_iBtstrp(:,1)<0 & RecordingSide_iBtstrp==1)|(i_selective & FR_pref_iBtstrp(:,1)>0 & RecordingSide_iBtstrp==2))/sum((i_selective & RecordingSide_iBtstrp==1)|(i_selective & RecordingSide_iBtstrp==2));
    
    % sample epoch
    i_selective = (sig_selective_iBtstrp(:,1));
    frac_contra_S_iBtstrp(i_btstrp,:) = sum((i_selective & FR_pref_iBtstrp(:,2)>0 & RecordingSide_iBtstrp==1)|(i_selective & FR_pref_iBtstrp(:,2)<0 & RecordingSide_iBtstrp==2))/sum((i_selective & RecordingSide_iBtstrp==1)|(i_selective & RecordingSide_iBtstrp==2));
    frac_ipsi_S_iBtstrp(i_btstrp,:) = sum((i_selective & FR_pref_iBtstrp(:,2)<0 & RecordingSide_iBtstrp==1)|(i_selective & FR_pref_iBtstrp(:,2)>0 & RecordingSide_iBtstrp==2))/sum((i_selective & RecordingSide_iBtstrp==1)|(i_selective & RecordingSide_iBtstrp==2));
    
    % delay epoch
    i_selective = (sig_selective_iBtstrp(:,2));
    frac_contra_D_iBtstrp(i_btstrp,:) = sum((i_selective & FR_pref_iBtstrp(:,3)>0 & RecordingSide_iBtstrp==1)|(i_selective & FR_pref_iBtstrp(:,3)<0 & RecordingSide_iBtstrp==2))/sum((i_selective & RecordingSide_iBtstrp==1)|(i_selective & RecordingSide_iBtstrp==2));
    frac_ipsi_D_iBtstrp(i_btstrp,:) = sum((i_selective & FR_pref_iBtstrp(:,3)<0 & RecordingSide_iBtstrp==1)|(i_selective & FR_pref_iBtstrp(:,3)>0 & RecordingSide_iBtstrp==2))/sum((i_selective & RecordingSide_iBtstrp==1)|(i_selective & RecordingSide_iBtstrp==2));
    
    % response epoch
    i_selective = (sig_selective_iBtstrp(:,3));
    frac_contra_R_iBtstrp(i_btstrp,:) = sum((i_selective & FR_pref_iBtstrp(:,4)>0 & RecordingSide_iBtstrp==1)|(i_selective & FR_pref_iBtstrp(:,4)<0 & RecordingSide_iBtstrp==2))/sum((i_selective & RecordingSide_iBtstrp==1)|(i_selective & RecordingSide_iBtstrp==2));
    frac_ipsi_R_iBtstrp(i_btstrp,:) = sum((i_selective & FR_pref_iBtstrp(:,4)<0 & RecordingSide_iBtstrp==1)|(i_selective & FR_pref_iBtstrp(:,4)>0 & RecordingSide_iBtstrp==2))/sum((i_selective & RecordingSide_iBtstrp==1)|(i_selective & RecordingSide_iBtstrp==2));
    
    
end

subplot(1,4,1); hold on
bar(1,frac_contra_SDR,'b');
bar(2,frac_ipsi_SDR,'r');
errorbar(1, frac_contra_SDR, std(frac_contra_SDR_iBtstrp),'b');
errorbar(2, frac_ipsi_SDR, std(frac_ipsi_SDR_iBtstrp),'r');
line([.5 2.5],[.5 .5],'color','k','linestyle',':')
ylim([0 1])
ylabel('frac of sig sel neurons')
legend('contra','ipsi');
title('whole trial')

subplot(1,4,2); hold on
bar(1,frac_contra_S,'b');
bar(2,frac_ipsi_S,'r');
errorbar(1, frac_contra_S, std(frac_contra_S_iBtstrp),'b');
errorbar(2, frac_ipsi_S, std(frac_ipsi_S_iBtstrp),'r');
line([.5 2.5],[.5 .5],'color','k','linestyle',':')
ylim([0 1])
title('sample')

subplot(1,4,3); hold on
bar(1,frac_contra_D,'b');
bar(2,frac_ipsi_D,'r');
errorbar(1, frac_contra_D, std(frac_contra_D_iBtstrp),'b');
errorbar(2, frac_ipsi_D, std(frac_ipsi_D_iBtstrp),'r');
line([.5 2.5],[.5 .5],'color','k','linestyle',':')
ylim([0 1])
title('delay')

subplot(1,4,4); hold on
bar(1,frac_contra_R,'b');
bar(2,frac_ipsi_R,'r');
errorbar(1, frac_contra_R, std(frac_contra_R_iBtstrp),'b');
errorbar(2, frac_ipsi_R, std(frac_ipsi_R_iBtstrp),'r');
line([.5 2.5],[.5 .5],'color','k','linestyle',':')
ylim([0 1])
title('response')




%% contra ipsi preference time course
time_step = [];
n_sig_contra = [];
n_sig_ipsi = [];
n_step = 0;
for i_timestep = -3:.1:2
    
    disp(['time step: ',num2str(i_timestep),' s']);
    
    n_step = n_step+1;
    
    n_sig_tmp = [];
    for i_unit = 1:size(spk_times_yes_correct_all,1)
        
        spike_times_yes = spk_times_yes_correct_all{i_unit,1};
        spike_times_no = spk_times_no_correct_all{i_unit,1};
        
        spike_count_yes = [];
        for i_trial = 1:size(spike_times_yes,1)
            spike_count_yes(i_trial,1) = sum(spike_times_yes{i_trial}>(i_timestep-.1) & spike_times_yes{i_trial}<i_timestep);
        end
        
        spike_count_no = [];
        for i_trial = 1:size(spike_times_no,1)
            spike_count_no(i_trial,1) = sum(spike_times_no{i_trial}>(i_timestep-.1) & spike_times_no{i_trial}<i_timestep);
        end
        
        [h p] = ttest2(spike_count_yes,spike_count_no);
        
        if p<.01
            if ((mean(spike_count_yes)>mean(spike_count_no)) & RecordingSide_all(i_unit,1)==1) | ((mean(spike_count_yes)<mean(spike_count_no)) & RecordingSide_all(i_unit,1)==2)
                n_sig_tmp(i_unit,1) = 1;    % contra
            elseif ((mean(spike_count_yes)<mean(spike_count_no)) & RecordingSide_all(i_unit,1)==1) | ((mean(spike_count_yes)>mean(spike_count_no)) & RecordingSide_all(i_unit,1)==2)
                n_sig_tmp(i_unit,1) = -1;    % ipsi
            else
                error('should not be here')
            end
        else
            n_sig_tmp(i_unit,1) = 0;
        end
        
    end
    
    time_step(n_step,1) = i_timestep;
    n_sig_contra(n_step,1) = sum(n_sig_tmp==1);
    n_sig_ipsi(n_step,1) = sum(n_sig_tmp==-1);
    
end

figure; hold on
bar(time_step,n_sig_contra,'b');
bar(time_step,-n_sig_ipsi,'r');
line([-2.6 -2.6],[-250 250],'color','k');
line([-1.3 -1.3],[-250 250],'color','k');
line([0 0],[-250 250],'color','k');
xlabel('Time (s)')
ylabel('Num of sig sel neurons')
legend('contra','ipsi')



