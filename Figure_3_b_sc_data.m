%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 3b # of selective SC neurons contra vs ipsi

clear all
close all
addpath('../func/');


load Figure_3_b_sc_data


sig_selective(isnan(sig_selective))=1;
sig_selective = sig_selective<0.01;%sig_selective<(0.05/size(sig_selective,1)/3); % <----------- This p value can be adjusted

cell_location = func_get_SC_core_shell(Unit_CCF(:,2), Unit_CCF(:,3), Unit_CCF(:,1));
i_SC_core = (cell_location==1);
i_SC_shell = (cell_location==2);


% ---------------- flip all units into the left hemisphere ------------
CCF_width = 1140; % size of CCF image, 10um volume
CCF_midline = 570;
i_units = find(Unit_CCF(:,2)>CCF_midline);
Unit_CCF(i_units,2) = CCF_width-Unit_CCF(i_units,2);        % Some sessions of rightSC recordings are already registered to left hemisphere in CCF. This is not reliable as measure to identify recorded hemisphere.

% Change label of 'yes' and 'no' trials, so 'yes' trials are also 'contra'
i_units = find(RecordingSide_all==2);   % flip 'yes' & 'no' for rightSC units, Ultimately, the direction is based on SC hemisphere indicated in sessionTagName

spk_count_no_tmp = spk_count_yes_screen(i_units,:);
spk_count_yes_tmp = spk_count_no_screen(i_units,:);
spk_count_yes_screen(i_units,:) = spk_count_yes_tmp;
spk_count_no_screen(i_units,:) = spk_count_no_tmp;


spk_count_no_tmp = spk_count_yes_all(i_units,:);
spk_count_yes_tmp = spk_count_no_all(i_units,:);
spk_count_yes_all(i_units,:) = spk_count_yes_tmp;
spk_count_no_all(i_units,:) = spk_count_no_tmp;


spk_count_no_tmp = spk_count_yes_error_all(i_units,:);
spk_count_yes_tmp = spk_count_no_error_all(i_units,:);
spk_count_yes_error_all(i_units,:) = spk_count_yes_tmp;
spk_count_no_error_all(i_units,:) = spk_count_no_tmp;


PSTH_no_tmp = PSTH_yes_cue_aligned(i_units,:);
PSTH_yes_tmp = PSTH_no_cue_aligned(i_units,:);
PSTH_yes_cue_aligned(i_units,:) = PSTH_yes_tmp;
PSTH_no_cue_aligned(i_units,:) = PSTH_no_tmp;

PSTH_no_tmp = PSTH_yes_error_cue_aligned(i_units,:);
PSTH_yes_tmp = PSTH_no_error_cue_aligned(i_units,:);
PSTH_yes_error_cue_aligned(i_units,:) = PSTH_yes_tmp;
PSTH_no_error_cue_aligned(i_units,:) = PSTH_no_tmp;


spike_times_tmp1 = spike_times_all(i_units,2);
spike_times_tmp2 = spike_times_all(i_units,1);
spike_times_tmp3 = spike_times_all(i_units,4);
spike_times_tmp4 = spike_times_all(i_units,3);

spike_times_all(i_units,1) = spike_times_tmp1;
spike_times_all(i_units,2) = spike_times_tmp2;
spike_times_all(i_units,3) = spike_times_tmp3;
spike_times_all(i_units,4) = spike_times_tmp4;



%% contra ipsi preference time course
n_sig_contra_all = [];
n_sig_ipsi_all = [];
for i_location = -1:2

    time_step = [];
    n_sig_contra = [];
    n_sig_ipsi = [];
    n_step = 0;
    for i_timestep = -3:.1:2
        
        disp(['time step: ',num2str(i_timestep),' s']);
        
        n_step = n_step+1;
        
        n_sig_tmp = [];
        n_unit = 0;
        for i_unit = find(cell_location == i_location)'%1:size(spike_times_all,1)
            
            n_unit = n_unit+1;
            
            spike_times_yes = spike_times_all{i_unit,1};
            spike_times_no = spike_times_all{i_unit,2};
            
            spike_count_yes = [];
            for i_trial = 1:size(spike_times_yes,1)
                spike_count_yes(i_trial,1) = sum(spike_times_yes{i_trial}>(i_timestep-.1) & spike_times_yes{i_trial}<(i_timestep+.1));
            end
            
            spike_count_no = [];
            for i_trial = 1:size(spike_times_no,1)
                spike_count_no(i_trial,1) = sum(spike_times_no{i_trial}>(i_timestep-.1) & spike_times_no{i_trial}<(i_timestep+.1));
            end
            
            [h p] = ttest2(spike_count_yes,spike_count_no);
            
            if p<.01
                if ((mean(spike_count_yes)>mean(spike_count_no)))
                    n_sig_tmp(n_unit,1) = 1;    % contra
                elseif ((mean(spike_count_yes)<mean(spike_count_no)))
                    n_sig_tmp(n_unit,1) = -1;    % ipsi
                else
                    error('should not be here')
                end
            else
                n_sig_tmp(n_unit,1) = 0;
            end
            
        end
        
        time_step(n_step,1) = i_timestep;
        n_sig_contra(n_step,1) = sum(n_sig_tmp==1);
        n_sig_ipsi(n_step,1) = sum(n_sig_tmp==-1);
        
    end
    
    n_sig_contra_all(end+1,:) = n_sig_contra;
    n_sig_ipsi_all(end+1,:) = n_sig_ipsi;

end

n_sig_contra_tmp = sum(n_sig_contra_all(2:4,:));
n_sig_ipsi_tmp = sum(n_sig_ipsi_all(2:4,:));

figure; hold on
bar(time_step,n_sig_contra_tmp,'b');
bar(time_step,-n_sig_ipsi_tmp,'r');
line([-2.6 -2.6],[-350 350],'color','k');
line([-1.3 -1.3],[-350 350],'color','k');
line([0 0],[-350 350],'color','k');
title('All SC')
ylim([-350 350])



