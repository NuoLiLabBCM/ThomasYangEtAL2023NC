%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 3b # of selective ALM neurons contra vs ipsi

clear all
close all
addpath('../func/');


load Figure_3_b_data


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
ylim([-450 450])
line([-2.6 -2.6],[-450 450],'LineStyle','--','color','k');
line([-1.3 -1.3],[-450 450],'LineStyle','--','color','k');
line([0 0],[-450 450],'LineStyle','--','color','k');

line([0 0],[-250 250],'color','k');
xlabel('Time (s)')
ylabel('Num of sig sel neurons')
legend('contra','ipsi')



