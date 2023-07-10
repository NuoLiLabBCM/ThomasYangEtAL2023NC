%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 7b selectivity of tagged GABAergic neurons

clear all
close all
addpath('../func')

load Figure_7_b_data

sig_selective = sig_selective_p<.01;        %<------- user defined


%% compute response to photostim and latency
spike_count_all = [];
latency_all = [];
jitter_all = [];
sig_response_all = [];
PSTH_tagging_all = [];youtu
spk_times_tagging_all = {};
PSTH_tagging_pulses_all = [];
waveform_mean_ctrl_all = [];
waveform_sem_ctrl_all = [];
waveform_mean_stim_all = [];
waveform_sem_stim_all = [];
for i_unit = 1:size(spk_times_stimPulse_all,1)
    
    spk_times_iUnit = spk_times_stimPulse_all{i_unit};
    stim_iUnit = stimPulse_condition_all{i_unit};
    waveform_iUnit = waveform_stimPulse_all{i_unit};
    
    % combine responses to all 3 pulses
    spk_times_tmp = {};
    spk_count_tmp = [];
    first_spk_tmp = [];
    for i_rep = 1:size(stim_iUnit,1)
        spk_times_tmp{end+1,1} = spk_times_iUnit{i_rep};
        spk_times_tmp{end+1,1} = spk_times_iUnit{i_rep}-.2;
        spk_times_tmp{end+1,1} = spk_times_iUnit{i_rep}-.4;
        
        spk_count_tmp(end+1,:) = [sum(spk_times_iUnit{i_rep}<0 & spk_times_iUnit{i_rep}>-.01)               sum(spk_times_iUnit{i_rep}>0 & spk_times_iUnit{i_rep}<.01)];
        spk_count_tmp(end+1,:) = [sum((spk_times_iUnit{i_rep}-.2)<0 & (spk_times_iUnit{i_rep}-.2)>-.01)     sum((spk_times_iUnit{i_rep}-.2)>0 & (spk_times_iUnit{i_rep}-.2)<.01)];
        spk_count_tmp(end+1,:) = [sum((spk_times_iUnit{i_rep}-.4)<0 & (spk_times_iUnit{i_rep}-.4)>-.01)     sum((spk_times_iUnit{i_rep}-.4)>0 & (spk_times_iUnit{i_rep}-.4)<.01)];
        
        i_1st_spk = find(spk_times_iUnit{i_rep}>0 & spk_times_iUnit{i_rep}<.01);
        if ~isempty(i_1st_spk)
            first_spk_tmp(end+1,:) = spk_times_iUnit{i_rep}(i_1st_spk(1));
        end
        i_1st_spk = find((spk_times_iUnit{i_rep}-.2)>0 & (spk_times_iUnit{i_rep}-.2)<.01);
        if ~isempty(i_1st_spk)
            first_spk_tmp(end+1,:) = spk_times_iUnit{i_rep}(i_1st_spk(1))-.2;
        end
        i_1st_spk = find((spk_times_iUnit{i_rep}-.2)>0 & (spk_times_iUnit{i_rep}-.2)<.01);
        if ~isempty(i_1st_spk)
            first_spk_tmp(end+1,:) = spk_times_iUnit{i_rep}(i_1st_spk(1))-.2;
        end
        
    end
    
    [psth1 t1] = func_getPSTH_smallBin(spk_times_tmp, -0.01, 0.05);%
    [h p] = ttest2(spk_count_tmp(:,1), spk_count_tmp(:,2));
    if p<0.01 & (mean(spk_count_tmp(:,2))>mean(spk_count_tmp(:,1)))
        [psth_peak i_latency] = max(psth1);
        t_latency = t1(i_latency(1));
    else
        t_latency = nan;
    end
    
    spike_count_all(i_unit,1) = mean(spk_count_tmp(:,2))-mean(spk_count_tmp(:,1));
    latency_all(i_unit,1) = t_latency;
    jitter_all(i_unit,:) = [mean(first_spk_tmp) std(first_spk_tmp)];
    sig_response_all(i_unit,1) = p;
    PSTH_tagging_all(i_unit,:) = psth1;
    spk_times_tagging_all{i_unit,1} = spk_times_tmp;
    
    
    % response to individal pulse
    [psth2 t2] = func_getPSTH_smallBin(spk_times_iUnit, -0.2, 0.8);
    PSTH_tagging_pulses_all(i_unit,:) = psth2;
    
    
    % waveform
    waveform_stim_tmp = [];
    waveform_ctrl_tmp = [];
    for i_rep = 1:size(stim_iUnit,1)
        i_pulse1 = (spk_times_iUnit{i_rep}>0 & spk_times_iUnit{i_rep}<.01);
        i_pulse2 = ((spk_times_iUnit{i_rep}-.2)>0 & (spk_times_iUnit{i_rep}-.2)<.01);
        i_pulse3 = ((spk_times_iUnit{i_rep}-.4)>0 & (spk_times_iUnit{i_rep}-.4)<.01);
        if ~isempty(i_pulse1 | i_pulse2 | i_pulse3)
            waveform_stim_tmp = cat(1,waveform_stim_tmp, waveform_iUnit{i_rep}(i_pulse1 | i_pulse2 | i_pulse3,:));
        end
        waveform_ctrl_tmp = cat(1,waveform_ctrl_tmp, waveform_iUnit{i_rep}(~(i_pulse1 | i_pulse2 | i_pulse3),:));
    end
    waveform_mean_ctrl_all(i_unit,:) = mean(waveform_ctrl_tmp);
    waveform_sem_ctrl_all(i_unit,:) = std(waveform_ctrl_tmp);%/sqrt(size(waveform_ctrl_tmp,1));
    waveform_mean_stim_all(i_unit,:) = mean(waveform_stim_tmp);
    waveform_sem_stim_all(i_unit,:) = std(waveform_stim_tmp);%/sqrt(size(waveform_stim_tmp,1));
    
end

figure; 
subplot(1,2,1); hold on
plot(latency_all*1000, spike_count_all,'o');
% set(gca,'xscale','log','yscale','log')
xlabel('latency (ms)')
ylabel('spike evoked per light pulse')

subplot(1,2,2); hold on
plot(latency_all*1000, jitter_all(:,2)*1000,'o');
% set(gca,'xscale','log','yscale','log')
xlabel('latency (ms)')
ylabel('jitter (ms)')

% for figures
i_discard = (isnan(latency_all));
figure
subplot(1,3,1);
[y x] = hist(spike_count_all(~i_discard),0:.2:2);
bar(x,y);
line([mean(spike_count_all(~i_discard)) mean(spike_count_all(~i_discard))],[0 max(y)*1.2])
xlabel('spike evoked per light pulse')
ylabel('# of neurons')

subplot(1,3,2);
[y x] = hist(latency_all(~i_discard)*1000,0:10);
bar(x,y);
line([mean(latency_all(~i_discard)*1000) mean(latency_all(~i_discard)*1000)],[0 max(y)*1.2])
xlabel('latency (ms)')
ylabel('# of neurons')

subplot(1,3,3);
plot(t2,mean(PSTH_tagging_pulses_all))
xlabel('time (s)')
ylabel('spk/s')
title('average response')

%% plot the PSTH of the tagged population
figure
subplot(1,4,1); hold on
i_pulse1 = find(t2>-.01 & t2<.05);
t_pulse1 = t2(i_pulse1);
bar(t_pulse1, mean(PSTH_tagging_pulses_all(:,i_pulse1)), 'k');
line([0 0],[0 1.2]*max(mean(PSTH_tagging_pulses_all)),'color','k')
xlim([-.01 .04])
title('Pulse 1')

subplot(1,4,2); hold on
i_pulse2 = find((t2-.2)>-.01 & (t2-.2)<.05);
t_pulse2 = t2(i_pulse2)-.2;
bar(t_pulse2, mean(PSTH_tagging_pulses_all(:,i_pulse2)), 'k');
line([0 0],[0 1.2]*max(mean(PSTH_tagging_pulses_all)),'color','k')
xlim([-.01 .04])
title('Pulse 2')

subplot(1,4,3); hold on
i_pulse3 = find((t2-.4)>-.01 & (t2-.4)<.05);
t_pulse3 = t2(i_pulse3)-.4;
bar(t_pulse3, mean(PSTH_tagging_pulses_all(:,i_pulse3)), 'k');
line([0 0],[0 1.2]*max(mean(PSTH_tagging_pulses_all)),'color','k')
xlim([-.01 .04])
title('Pulse 3')


subplot(1,4,4); hold on
bar(t1, mean(PSTH_tagging_all), 'k');
line([0 0],[0 1.2]*max(mean(PSTH_tagging_all)),'color','k')
xlim([-.01 .04])
title('All pulses combined')



%% plot individual tagged neurons
% plot the selectivity and PSTH of the tagged units
figure; subplot(4,6,1); title('tagged');
n_plot = 0;
for i_unit = 1:size(PSTH_yes_correct_cue_aligned,1)
    
    n_plot = n_plot+1;
    if n_plot>24
        n_plot=1;
        figure; subplot(4,6,1); title('tagged');
    end
    
    if RecordingSide_all(i_unit,1) == 1
        psth1 = PSTH_yes_correct_cue_aligned(i_unit,:);
        psth2 = PSTH_no_correct_cue_aligned(i_unit,:);
        spike_times1 = spike_times_all{i_unit,1};       % yes correct trials
        spike_times2 = spike_times_all{i_unit,2};       % no correct trials  
    elseif RecordingSide_all(i_unit,1) == 2
        psth2 = PSTH_yes_correct_cue_aligned(i_unit,:);
        psth1 = PSTH_no_correct_cue_aligned(i_unit,:);
        spike_times2 = spike_times_all{i_unit,2};       % no correct trials  
        spike_times1 = spike_times_all{i_unit,1};       % yes correct trials
    else
        error('d');
    end
    
    subplot(4,6,n_plot); hold on
    plot(t,psth2,'r');
    plot(t,psth1,'b');
    line([0 0],[0 2.2]*max([psth1 psth2]),'color','k');
    line([-1.3 -1.3],[0 2.2]*max([psth1 psth2]),'color','k');
    line([-2.6 -2.6],[0 2.2]*max([psth1 psth2]),'color','k');
    xlim([-3 2]);
    ylim([0 2.3]*max([psth1 psth2]));
    
    filename_tmp = Filename_all{i_unit};
    istr = findstr(filename_tmp,'_allData');
    title(filename_tmp(1:istr(1)));
    
end



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
[h p_conta_sel_delay] = ttest(mean(contra_sel(:,i_delay),2),0,'tail','left');

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
for i_btstrp = 1:10000
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
title(['Bootstrp across mice (n=',num2str(length(unique(Mice_all(:,1)))),'); p=',num2str(sum(frac_contra_ipsi_btstrp(:,1)>frac_contra_ipsi_btstrp(:,2))/size(frac_contra_ipsi_btstrp,1))]);




