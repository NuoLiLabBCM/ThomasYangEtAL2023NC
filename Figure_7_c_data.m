%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 7c selectivity of tagged glutamatergic neurons

clear all
close all
addpath('../func')

load Figure_7_c_data

%% plot contra selectivity
i_sel_unit = (sig_selective(:,1)==1 | sig_selective(:,2)==1);

i_t = find(t(1,:)<-.1);

figure; hold on
psth1 = PSTH_yes_correct_cue_aligned(i_sel_unit,:);
psth2 = PSTH_no_correct_cue_aligned(i_sel_unit,:);
side = RecordingSide_all(i_sel_unit,:);
contra_sel = psth1-psth2;
contra_sel(side==2,:) = -contra_sel(side==2,:);
func_plot_mean_and_sem(t(i_t),contra_sel(:,i_t), 'k', [.6 .6 .6], 'n');

i_delay = find(t(1,:)>-1.4 & t(1,:)<-.1);
[h p_conta_sel_delay] = ttest(mean(contra_sel(:,i_delay),2),0,'tail','right');

line([-3 1.2],[0 0],'color','k')
line([0 0],[-8 8],'color','k')
line([-1.3 -1.3],[-8 8],'color','k')
line([-2.6 -2.6],[-8 8],'color','k')
xlabel('Time')
ylabel('Contra selectivity (spk/s)')
legend('tagged','tagged','non-tagged','non-tagged')
title(['Tagged and selective (n=',num2str(sum(i_sel_unit)),'), delay ipsi sel, p=',num2str(p_conta_sel_delay)]);




% spike counts used to determine the preferred trial type
spk_count_yes_all = [];
spk_count_no_all = [];
for i_unit = 1:size(spike_times_all,1)
   
    spike_times_yes = spike_times_all{i_unit,1};
    spike_times_no = spike_times_all{i_unit,2};

    % yes trials correct
    spk_count_yes = [];
    n_trial = 0;
    for i_trial = 1:size(spike_times_yes,1)
        spk_count_tmp = [
            sum(spike_times_yes{i_trial,1}>-2.6 & spike_times_yes{i_trial,1}<1.3),... % whole trial
            sum(spike_times_yes{i_trial,1}>-2.6 & spike_times_yes{i_trial,1}<-1.3),... % sample
            sum(spike_times_yes{i_trial,1}>-1.3 & spike_times_yes{i_trial,1}<0),... % delay
            sum(spike_times_yes{i_trial,1}>0 & spike_times_yes{i_trial,1}<1.3),... % response
            ];
        spk_count_yes(end+1,:) = spk_count_tmp;
    end
    clear spk_count_tmp

    
    % no trials  screen
    spk_count_no = [];
    n_trial = 0;
    for i_trial = 1:size(spike_times_no,1)
        spk_count_tmp = [
            sum(spike_times_no{i_trial,1}>-2.6 & spike_times_no{i_trial,1}<1.3),... % whole trial
            sum(spike_times_no{i_trial,1}>-2.6 & spike_times_no{i_trial,1}<-1.3),... % sample
            sum(spike_times_no{i_trial,1}>-1.3 & spike_times_no{i_trial,1}<0),... % delay
            sum(spike_times_no{i_trial,1}>0 & spike_times_no{i_trial,1}<1.3),... % response
            ];
        spk_count_no(end+1,:) = spk_count_tmp;
    end
    clear spk_count_tmp

    spk_count_yes_all(i_unit,:) = mean(spk_count_yes);
    spk_count_no_all(i_unit,:) = mean(spk_count_no);

end


% fraction contra vs. ipsi
figure; hold on
FR_pref = spk_count_yes_all-spk_count_no_all;
FR_pref(RecordingSide_all==2,:) = -FR_pref(RecordingSide_all==2,:);

i_sel_unit = (sig_selective(:,2)==1 | sig_selective(:,1)==1 );
bar(1,sum(FR_pref(i_sel_unit,3)>0)/sum(i_sel_unit),'b')
bar(2,sum(FR_pref(i_sel_unit,3)<0)/sum(i_sel_unit),'r')


% bootstrp
frac_contra_ipsi_btstrp = [];
for i_btstrp = 1:100
    if rem(i_btstrp,1000)==0
        i_btstrp
    end
    
    Mice_ID = unique(Mice_all(:,1));
    Mice_ID_iBtstrp = randsample(Mice_ID,length(Mice_ID),'true');

    FR_pref_iBtstrp = [];
    sig_selective_iBtstrp = [];
    for i_mice = Mice_ID_iBtstrp'
        FR_pref_iBtstrp = [FR_pref_iBtstrp; FR_pref(Mice_all(:,1)==i_mice,:)];
        sig_selective_iBtstrp = [sig_selective_iBtstrp; sig_selective(Mice_all(:,1)==i_mice,:)];
    end
    
    i_sel_unit_iBtstrp = (sig_selective_iBtstrp(:,2)==1 | sig_selective_iBtstrp(:,1)==1 );
    frac_contra_ipsi_btstrp(i_btstrp,:) = [sum(FR_pref_iBtstrp(i_sel_unit_iBtstrp,3)>0)/sum(i_sel_unit_iBtstrp)   sum(FR_pref_iBtstrp(i_sel_unit_iBtstrp,3)<0)/sum(i_sel_unit_iBtstrp)];
end

errorbar(1,sum(FR_pref(i_sel_unit,3)>0)/sum(i_sel_unit),std(frac_contra_ipsi_btstrp(:,1)),'k');
scatter(ones(size(frac_contra_ipsi_btstrp,1),1),frac_contra_ipsi_btstrp(:,1),36,[0.9 0.9 0.9]);
errorbar(2,sum(FR_pref(i_sel_unit,3)<0)/sum(i_sel_unit),std(frac_contra_ipsi_btstrp(:,2)),'k');
scatter(2*ones(size(frac_contra_ipsi_btstrp,1),1),frac_contra_ipsi_btstrp(:,2),36,[0.9 0.9 0.9]);

ylim([0 1])
yline(0.5,':');
title(['Bootstrp across mice (n=',num2str(length(unique(Mice_all(:,1)))),'); p=',num2str(sum(frac_contra_ipsi_btstrp(:,1)<frac_contra_ipsi_btstrp(:,2))/size(frac_contra_ipsi_btstrp,1))]);





%%


sig_selective = sig_selective_p<.01;        %<------- user defined


% ---------------- flip all units into the left hemisphere ------------
% Change label of 'yes' and 'no' trials, so 'yes' trials are also 'contra'
i_units = find(RecordingSide_all==2);   % flip 'yes' & 'no' for rightSC units, Ultimately, the direction is based on SC hemisphere indicated in sessionTagName

PSTH_no_tmp = PSTH_yes_correct_cue_aligned(i_units,:);
PSTH_yes_tmp = PSTH_no_correct_cue_aligned(i_units,:);
PSTH_yes_correct_cue_aligned(i_units,:) = PSTH_yes_tmp;
PSTH_no_correct_cue_aligned(i_units,:) = PSTH_no_tmp;

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



%% contra vs. ipsi

    
time_step = [];
n_sig_contra = [];
n_sig_ipsi = [];
n_step = 0;
for i_timestep = -3:.1:2
    
    disp(['time step: ',num2str(i_timestep),' s']);
    
    n_step = n_step+1;
    
    n_sig_tmp = [];
    n_unit = 0;
    for i_unit = 1:size(spike_times_all,1)
        
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


figure; hold on
line([0 0],[-20 20],'color','k')
line([0 0]-2.6,[-20 20],'color','k')
line([0 0]-1.3,[-20 20],'color','k')
bar(time_step,n_sig_contra,'b');
bar(time_step,-n_sig_ipsi,'r');
xlim([-3 2])

xlabel('time of 1st lick (s)')
ylabel('# of Neurons')




%% plot contra selectivity
i_sel_unit = (sig_selective(:,1)==1 | sig_selective(:,2)==1);

i_t = find(t(1,:)<2);

figure; hold on
psth1 = PSTH_yes_correct_cue_aligned(i_sel_unit,:);
psth2 = PSTH_no_correct_cue_aligned(i_sel_unit,:);
side = RecordingSide_all(i_sel_unit,:);
contra_sel = psth1-psth2;
func_plot_mean_and_sem(t(i_t),contra_sel(:,i_t), 'k', [.6 .6 .6], 'n');

i_delay = find(t(1,:)>-1.4 & t(1,:)<-.1);
[h p_conta_sel_delay] = ttest(mean(contra_sel(:,i_delay),2),0,'tail','left');

line([-3 1.2],[0 0],'color','k')
line([0 0],[-8 8],'color','k')
line([-1.3 -1.3],[-8 8],'color','k')
line([-2.6 -2.6],[-8 8],'color','k')
xlim([-3 2])
xlabel('Time')
ylabel('Contra selectivity (spk/s)')
legend('tagged','tagged','non-tagged','non-tagged')
title(['Tagged and selective (n=',num2str(sum(i_sel_unit)),'), delay ipsi sel, p=',num2str(p_conta_sel_delay)]);



