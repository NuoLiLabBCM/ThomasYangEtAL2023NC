%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 2e ALM selectivity map

clear all
close all
addpath('../func/');

load Figure_2_e_data


im_CCF = imread(['..\CCF_alignment\ARA_template\AllenRefVolCoronal_100275.tif']);

im_CCF_anatomy = imread(['.\ALM2x_mergedRGB_slice300_cltr1247_P55_NL09.tif']);
im_CCF_anatomy = imresize(im_CCF_anatomy,2);


%% plot contra vs. ipsi selectivity 
sig_selective(isnan(sig_selective))=1;
% sig_selective = sig_selective<.01; % <----------- This p value can be adjusted
sig_selective = sig_selective<(0.05/size(sig_selective,1)/3); % <----------- This p value can be adjusted


selectivity_all = [];
for i_unit = 1:size(CCF_x_y_slice,1)
    
    
    spk_times_yes_tmp  = spk_times_yes_correct_all{i_unit};
    spk_times_no_tmp  = spk_times_no_correct_all{i_unit};
    
    Delay_dur_tmp = Delay_dur_all(i_unit,1);
    
    spk_count_yes = [];
    for i_rep = 1:size(spk_times_yes_tmp,1)
        spk_count_yes(i_rep,1) = sum(spk_times_yes_tmp{i_rep}<-Delay_dur_tmp & spk_times_yes_tmp{i_rep}>(-Delay_dur_tmp-1.3))/1.3;
        spk_count_yes(i_rep,2) = sum(spk_times_yes_tmp{i_rep}<0 & spk_times_yes_tmp{i_rep}>(-Delay_dur_tmp))/Delay_dur_tmp;
        spk_count_yes(i_rep,3) = sum(spk_times_yes_tmp{i_rep}<1.3 & spk_times_yes_tmp{i_rep}>0)/1.3;
    end
    
    spk_count_no = [];
    for i_rep = 1:size(spk_times_no_tmp,1)
        spk_count_no(i_rep,1) = sum(spk_times_no_tmp{i_rep}<-Delay_dur_tmp & spk_times_no_tmp{i_rep}>(-Delay_dur_tmp-1.3))/1.3;
        spk_count_no(i_rep,2) = sum(spk_times_no_tmp{i_rep}<0 & spk_times_no_tmp{i_rep}>(-Delay_dur_tmp))/Delay_dur_tmp;
        spk_count_no(i_rep,3) = sum(spk_times_no_tmp{i_rep}<1.3 & spk_times_no_tmp{i_rep}>0)/1.3;
    end
    
    selectivity_all(i_unit,:) = mean(spk_count_yes)-mean(spk_count_no);
end

% selectivity_all(selectivity_all>5)=5;

figure;
for i_tmp = 1:3
    subplot(1,3,i_tmp);
    image(im_CCF_anatomy); colormap gray; hold on
    plot(CCF_x_y_slice(:,1),CCF_x_y_slice(:,2),'ow','markersize',1.5)
    xlim([250 570])
    ylim([150 350]);
end

for i_unit = 1:size(CCF_x_y_slice,1)
    
    % sample
    my_dot_size = abs(selectivity_all(i_unit,1))/20*8+2;     % selectity amplitude
    my_dot_size(my_dot_size>20)=10;             % cap at 20
    if sig_selective(i_unit,1)>0
        subplot(1,3,1); hold on
        if selectivity_all(i_unit,1)>0
            plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'ow','markerfacecolor','b','markersize',my_dot_size);
        elseif selectivity_all(i_unit,1)<0
            plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'ow','markerfacecolor','r','markersize',my_dot_size);
        end
    end
    
    %delay
    my_dot_size = abs(selectivity_all(i_unit,2))/20*8+2;     % selectity amplitude
    my_dot_size(my_dot_size>20)=10;             % cap at 20
    if sig_selective(i_unit,2)>0
        subplot(1,3,2);
        if selectivity_all(i_unit,2)>0
            plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'ow','markerfacecolor','b','markersize',my_dot_size);
        elseif selectivity_all(i_unit,2)<0
            plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'ow','markerfacecolor','r','markersize',my_dot_size);
        end
    end
    
    %response
    my_dot_size = abs(selectivity_all(i_unit,3))/20*8+2;     % selectity amplitude
    my_dot_size(my_dot_size>20)=10;             % cap at 20
    if sig_selective(i_unit,3)>0
        subplot(1,3,3);
        if selectivity_all(i_unit,3)>0
            plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'ow','markerfacecolor','b','markersize',my_dot_size);
        elseif selectivity_all(i_unit,3)<0
            plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'ow','markerfacecolor','r','markersize',my_dot_size);
        end
    end
    
end



%% stimulus, choice, movement selectivity, by epoch
p_coeff_alpha = []; 
p_coeff_beta = []; 
p_coeff_sigma = []; 
for i_unit = 1:size(CCF_x_y_slice,1)
    
    spike_times1 = spk_times_yes_correct_all{i_unit,1};       %[yes no yes_error no_error]
    spike_times2 = spk_times_no_correct_all{i_unit,1};       %[yes no yes_error no_error]
    spike_times3 = spk_times_yes_error_all{i_unit,1};       %[yes no yes_error no_error]
    spike_times4= spk_times_no_error_all{i_unit,1};       %[yes no yes_error no_error]

    alpha = p_coeff_all(i_unit,:,1);
    beta = p_coeff_all(i_unit,:,2);
    sigma = p_coeff_all(i_unit,:,3);
    
    if ~isempty(spike_times3) & ~isempty(spike_times4)       %[stimulus, delay, response, whole trials, sample+delay]
        
        p_coeff_alpha(i_unit,1:5) = alpha;
        p_coeff_beta(i_unit,1:5) = beta;
        p_coeff_sigma(i_unit,1:5) = sigma;

    else
        p_coeff_alpha(i_unit,1:5) = [nan nan nan nan nan];
        p_coeff_beta(i_unit,1:5) = [nan nan nan nan nan];
        p_coeff_sigma(i_unit,1:5) = [nan nan nan nan nan];
        
    end    
    
end
    

max_sel = max(max(abs(-log10([p_coeff_alpha; p_coeff_beta; p_coeff_sigma]))));


% plot on CCF
% stimulus selectivity
figure;
for i_tmp = 1:3
    subplot(3,3,i_tmp);
    image(im_CCF); colormap gray; hold on
    plot(CCF_x_y_slice(:,1),CCF_x_y_slice(:,2),'ow','markersize',1)
    xlim([250 570])
    ylim([150 350]);
end
for i_unit = 1:size(CCF_x_y_slice,1)
    
    % sample
    if p_coeff_alpha(i_unit,1)<0.01
        subplot(3,3,1);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'og','markerfacecolor','g','markersize', -log10(p_coeff_alpha(i_unit,1))/max_sel*20+2);
    end
    
    %delay
    if p_coeff_alpha(i_unit,2)<0.01
        subplot(3,3,2);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'og','markerfacecolor','g','markersize', -log10(p_coeff_alpha(i_unit,2))/max_sel*20+2);
    end
    
    %response
    if p_coeff_alpha(i_unit,3)<0.01
        subplot(3,3,3);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'og','markerfacecolor','g','markersize', -log10(p_coeff_alpha(i_unit,3))/max_sel*20+2);
    end
    
end
subplot(3,3,1); title('Stimulus, Sample')
subplot(3,3,2); title('Delay')
subplot(3,3,3); title('Response')



% choice selectivity
for i_tmp = 1:3
    subplot(3,3,i_tmp+3);
    image(im_CCF); colormap gray; hold on
    plot(CCF_x_y_slice(:,1),CCF_x_y_slice(:,2),'ow','markersize',1)
    xlim([250 570])
    ylim([150 350]);
end
for i_unit = 1:size(CCF_x_y_slice,1)
    
    % sample
    if p_coeff_beta(i_unit,1)<0.01
        subplot(3,3,1+3);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'om','markerfacecolor','m','markersize', -log10(p_coeff_beta(i_unit,1))/max_sel*20+2);
    end
    
    %delay
    if p_coeff_beta(i_unit,2)<0.01
        subplot(3,3,2+3);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'om','markerfacecolor','m','markersize', -log10(p_coeff_beta(i_unit,2))/max_sel*20+2);
    end
    
    %response
    if p_coeff_beta(i_unit,3)<0.01
        subplot(3,3,3+3);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'om','markerfacecolor','m','markersize', -log10(p_coeff_beta(i_unit,3))/max_sel*20+2);
    end
    
end
subplot(3,3,4); title('Choice')



% outcome selectivity
for i_tmp = 1:3
    subplot(3,3,i_tmp+6);
    image(im_CCF); colormap gray; hold on
    plot(CCF_x_y_slice(:,1),CCF_x_y_slice(:,2),'ow','markersize',1)
    xlim([250 570])
    ylim([150 350]);
end

for i_unit = 1:size(CCF_x_y_slice,1)
    
    % sample
    if p_coeff_sigma(i_unit,1)<0.01
        subplot(3,3,1+6);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'oc','markerfacecolor','c','markersize', -log10(p_coeff_sigma(i_unit,1))/max_sel*20+2);
    end
    
    %delay
    if p_coeff_sigma(i_unit,2)<0.01
        subplot(3,3,2+6);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'oc','markerfacecolor','c','markersize', -log10(p_coeff_sigma(i_unit,2))/max_sel*20+2);
    end
    
    %response
    if p_coeff_sigma(i_unit,3)<0.01
        subplot(3,3,3+6);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'oc','markerfacecolor','c','markersize', -log10(p_coeff_sigma(i_unit,3))/max_sel*20+2);
    end
    
end
subplot(3,3,4); title('Outcome')



%% Licking Movement activity
lick_sel_cell = Lick_spike_count(:,6)<(.01/size(Lick_spike_count,1));

Mod_all = abs(Lick_spike_count(:,5));

% Licking modulation
figure;
image(im_CCF); colormap gray; hold on
plot(CCF_x_y_slice(:,1),CCF_x_y_slice(:,2),'ow','markersize',1)
xlim([250 570])
ylim([150 350]);

for i_unit = 1:size(CCF_x_y_slice,1)
    
    if lick_sel_cell(i_unit,1)
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'oy','markerfacecolor','y','markersize', Mod_all(i_unit)*20);
        %plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'oy','markerfacecolor','y','markersize', -log10(Lick_spike_count(i_unit,6))/max_sel*20+.1);
    end
    
end
title('Response epoch Licking modulation')






%% Detailed map of Stimulus, Choice, Movement selectivity MAP (for these maps, the most posterior neurons are exclude because they are outside of ALM and including them make them appear to be outside the brain).
CCF_sections_all = 250:20:320;
i_unit_section = [];
for i_unit = 1:size(CCF_x_y_slice,1)
    [dummy i_CCF] = min(abs(CCF_sections_all-CCF_x_y_slice(i_unit,3)));
    i_unit_section(i_unit,1)=i_CCF;
end

figure;
for i_section = 1:length(CCF_sections_all)    
    
    i_sel_unit = find(i_unit_section==i_section);
    
    % CCF image
    im_CCF = imread(['..\CCF_alignment\ARA_template\AllenRefVolCoronal_100',num2str(CCF_sections_all(i_section)),'.tif']);
    
    subplot(3,4,i_section);
    imagesc(im_CCF); hold on
    colormap('gray')
    plot(CCF_x_y_slice(i_sel_unit,1),CCF_x_y_slice(i_sel_unit,2),'ow','markersize',1)
    xlim([300 520]);    ylim([150 350]);
    title(['Section ',num2str(CCF_sections_all(i_section))])
    
    
    % stimulus selectivity
    %max_sel = max(max(abs(-log10([p_coeff_alpha(:,5)]))));
    selectivity_stim = abs(spike_count_yes_stim - spike_count_no_stim);
    for i_unit = i_sel_unit'
        
        % sample + delay
        if p_coeff_alpha(i_unit,5)<0.01
            subplot(3,4,i_section);
            %plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'og','markerfacecolor','g','markersize', -log10(p_coeff_alpha(i_unit,5))/max_sel*20);
            plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'og','markerfacecolor','g','markersize', selectivity_stim(i_unit,5));
        end
        
    end
    
    
    
    subplot(3,4,i_section+4);
    imagesc(im_CCF); hold on
    colormap('gray')
    plot(CCF_x_y_slice(i_sel_unit,1),CCF_x_y_slice(i_sel_unit,2),'ow','markersize',1)
    xlim([300 520]);    ylim([150 350]);
    
    % choice selectivity
    max_sel = max(max(abs(-log10([p_coeff_beta(:,5)]))));
    selectivity_choice = abs(spike_count_lickYes_choice - spike_count_lickNo_choice);
    for i_unit = i_sel_unit'
        
        % sample + delay
        if p_coeff_beta(i_unit,5)<0.01
            subplot(3,4,i_section+4);
            %plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'om','markerfacecolor','m','markersize', -log10(p_coeff_beta(i_unit,5))/max_sel*20);
            plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'om','markerfacecolor','m','markersize', selectivity_choice(i_unit,5));
        end
        
    end

    

    subplot(3,4,i_section+8);
    imagesc(im_CCF); hold on
    colormap('gray')
    plot(CCF_x_y_slice(i_sel_unit,1),CCF_x_y_slice(i_sel_unit,2),'ow','markersize',1)
    xlim([300 520]);    ylim([150 350]);
    
    % Licking modulation
    %max_sel = max(max(abs(-log10([Lick_spike_count(:,6)/size(Lick_spike_count,1)]))));
    for i_unit = i_sel_unit'
        
        if Lick_spike_count(i_unit,6)<(.01/size(Lick_spike_count,1))
            subplot(3,4,i_section+8);
            %plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'oy','markerfacecolor','y','markersize', -log10(Lick_spike_count(i_unit,6)/size(Lick_spike_count,1))/max_sel*20);
            plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'oy','markerfacecolor','y','markersize', Mod_all(i_unit)*10);
        end
        
    end

    
end





%% Summary Stimulus, Choice, Movement selectivity MAP (for these maps, the most posterior neurons are exclude because they are outside of ALM and including them make them appear to be outside the brain).
im_CCF = imread(['..\CCF_alignment\ARA_template\AllenRefVolCoronal_100275.tif']);

figure;
for i_tmp = 1:3
    subplot(3,1,i_tmp);
    image(im_CCF); colormap gray; hold on
    plot(CCF_x_y_slice(CCF_x_y_slice(:,3)<320,1),CCF_x_y_slice(CCF_x_y_slice(:,3)<320,2),'ow','markersize',1)
    xlim([250 570])
    ylim([150 350]);
end

% stimulus selectivity
selectivity_stim = abs(spike_count_yes_stim - spike_count_no_stim);
for i_unit = 1:size(CCF_x_y_slice,1)
    
    % sample + delay
    if p_coeff_alpha(i_unit,5)<0.01 & CCF_x_y_slice(i_unit,3)<320
        subplot(3,1,1);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'og','markerfacecolor','g','markersize', selectivity_stim(i_unit,5));
    end
    
end
subplot(3,1,1); title('Stimulus selectivity, Sample + Delay')

% choice selectivity
max_sel = max(max(abs(-log10([p_coeff_beta(:,5)]))));
selectivity_choice = abs(spike_count_lickYes_choice - spike_count_lickNo_choice);
for i_unit = 1:size(CCF_x_y_slice,1)
    
    % sample + delay
    if p_coeff_beta(i_unit,5)<0.01 & CCF_x_y_slice(i_unit,3)<320
        subplot(3,1,2);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'om','markerfacecolor','m','markersize', selectivity_choice(i_unit,5));
    end
    
end
subplot(3,1,2); title('Choice selectivity, Sample + Delay')


% Licking modulation
for i_unit = 1:size(Lick_spike_count,1)
    
    if Lick_spike_count(i_unit,6)<(.01/size(Lick_spike_count,1)) & CCF_x_y_slice(i_unit,3)<320
        subplot(3,1,3);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'oy','markerfacecolor','y','markersize', Mod_all(i_unit)*20);
    end
    
end
title('Licking modulation, Response')




% ---------- plot selectivity time course -----------------
% align stimulus selectivity by sample epoch start
PSTH_yes_stim_aligned = [];
PSTH_no_stim_aligned = [];
for i_cell = 1:size(PSTH_yes_stim,1)  
    disp(i_cell)
    i_t = t>(-Delay_dur_all(i_cell,1)-1.3-.3) & t<(4.5-Delay_dur_all(i_cell,1)-1.3-.3);
    t_stim = t(i_t)+Delay_dur_all(i_cell,1)+1.3;
    PSTH_yes_stim_aligned(i_cell,:) = PSTH_yes_stim(i_cell,i_t);
    PSTH_no_stim_aligned(i_cell,:) = PSTH_no_stim(i_cell,i_t);
end


% stimuluus selectivity
selectivity_stim = PSTH_yes_stim_aligned - PSTH_no_stim_aligned;  
i_sel_cell = spike_count_yes_stim(:,1)<spike_count_no_stim(:,1);
selectivity_stim(i_sel_cell,:) = -selectivity_stim(i_sel_cell,:);

figure;
i_sel_cell = ~isnan(sum(selectivity_stim,2)) & p_coeff_alpha(:,1)<0.01;
subplot(3,1,1); hold on
func_plot_mean_and_sem(t_stim, selectivity_stim(i_sel_cell,:), 'g', [.7 1 .7], 'n');
ylim([-.5 5])
line([0 0],[-.5 5],'color','k')
xlim([-1 5]);

% choice selectivity
selectivity_choice = PSTH_lickYes_choice - PSTH_lickNo_choice; 
i_sel_cell = spike_count_lickYes_choice(:,5)<spike_count_lickNo_choice(:,5);
selectivity_choice(i_sel_cell,:) = -selectivity_choice(i_sel_cell,:);

i_sel_cell = ~isnan(sum(selectivity_choice,2)) & p_coeff_beta(:,5)<0.01;
subplot(3,1,2); hold on
func_plot_mean_and_sem(t, selectivity_choice(i_sel_cell,:), 'm', [1 .7 1], 'n');
ylim([-.5 6])
line([0 0],[-.5 6],'color','k')



% licking activity
i_sel_cell = find(Lick_spike_count(:,6)<(.01/size(Lick_spike_count,1)));

selectivity_lick = [];
activity_modulation_lick = [];
n_cell = 0;
for i_cell = i_sel_cell'
    
    n_cell = n_cell+1;
    
    % all trials, correct
    spike_times_psth = {};
    spike_times_tmp = [spk_times_yes_correct_all{i_cell,1}; spk_times_no_correct_all{i_cell,1}];
    lick_times_tmp = [lick_times_all{i_cell,1}; lick_times_all{i_cell,2}];
    n_trial = 0;
    for i_trial = 1:size(spike_times_tmp,1)
        n_trial = n_trial+1;
        time_1stLick = lick_times_tmp{i_trial};
        time_1stLick = time_1stLick(time_1stLick(:,2)>0,2);
        time_1stLick = time_1stLick(1);
        spike_times_psth{n_trial,1} = spike_times_tmp{i_trial} - time_1stLick;
    end
    [psth t_lick1] = func_getPSTH_smallBin(spike_times_psth,-.5,1);
    
    selectivity_lick(n_cell,:) = psth/max(psth);
    
    
    % activity modulation by licking
    spike_times_tmp = [spk_times_yes_correct_all{i_cell,1}; spk_times_no_correct_all{i_cell,1}];
    lick_times_tmp = [lick_times_all{i_cell,1}; lick_times_all{i_cell,2}];
    
    n_trials_tmp = length(spike_times_tmp);
    i_lick_trial = randsample(n_trials_tmp,n_trials_tmp);
    i_lick_trial_screen = i_lick_trial(1:20);
    i_lick_trial = i_lick_trial(21:end);
    
    spike_times_psth = {};
    spike_count_lick = [];      %[baseline   sample   delay   response]
    n_trial = 0;
    for i_trial = 1:size(spike_times_tmp,1)
        n_trial = n_trial+1;
        time_1stLick = lick_times_tmp{i_trial};
        time_1stLick = time_1stLick(time_1stLick(:,2)>0,2);
        time_1stLick = time_1stLick(1);
        spike_times_psth{n_trial,1} = spike_times_tmp{i_trial} - time_1stLick;
        
        Delay_dur_tmp = Delay_dur_all(i_cell,1);
        
        spike_count_tmp(n_trial,:) = [
            sum(spike_times_tmp{i_trial}>(-Delay_dur_tmp-1.3-.5) & spike_times_tmp{i_trial}<(-Delay_dur_tmp-1.3)),...           % baseline
            sum(spike_times_tmp{i_trial}>(-Delay_dur_tmp-1.3) & spike_times_tmp{i_trial}<-Delay_dur_tmp),...                    % sample
            sum(spike_times_tmp{i_trial}>-Delay_dur_tmp & spike_times_tmp{i_trial}<0),...     % delay
            sum(spike_times_tmp{i_trial}>0 & spike_times_tmp{i_trial}<1.3),...  % response
            ];
        
    end
    [psth t_lick2] = func_getPSTH_smallBin(spike_times_psth(i_lick_trial,1),-1.5,1.5);
    spike_count_lick_screen = mean(spike_count_tmp(i_lick_trial_screen,:));
    
    if spike_count_lick_screen(4)<(spike_count_lick_screen(3)+spike_count_lick_screen(2))
        activity_modulation_lick(n_cell,:) = -psth+mean(psth(t_lick2<-.5));
    else
        activity_modulation_lick(n_cell,:) = psth-mean(psth(t_lick2<-.5));
    end
    
end
subplot(3,1,3); hold on
line([0 0],[-2 6],'color','k')
func_plot_mean_and_sem(t_lick2(t_lick2>-1.2), activity_modulation_lick(:,t_lick2>-1.2), [.8 .8 0], [1 1 0], 'n');
ylim([-2 6])
xlim([-1.5 0.8])

figure; %hold on
imagesc(t_lick1, linspace(1,1,size(selectivity_lick,2)), selectivity_lick);
line([0 0],[0 208],'color','y','linewidth',1)
% colormap pink



%% plot activity of licking modulated neurons
i_sel_cell = find(Lick_spike_count(:,6)<(.01/size(Lick_spike_count,1)));

Mod_all = abs(Lick_spike_count(:,1))+abs(Lick_spike_count(:,3));
[dummy, i_sort] = sort(Mod_all(i_sel_cell),'descend');
i_sel_cell = i_sel_cell(i_sort);

figure
n_cell = 0;
for i_cell = i_sel_cell(1:30)'      %<---------- plot top 30 cells only
    
    n_cell = n_cell+1;
    
    if n_cell>30
        figure;
        n_cell = 1;
    end
    
    % yes trials, correct
    spike_times_psth = {};
    spk_times_tmp = spk_times_yes_correct_all{i_cell,1};
    lick_times_tmp = lick_times_all{i_cell,1};
    n_trial = 0;
    for i_trial = 1:size(spk_times_tmp,1)
        n_trial = n_trial+1;
        time_1stLick = lick_times_tmp{i_trial};
        time_1stLick = time_1stLick(time_1stLick(:,2)>0,2);
        time_1stLick = time_1stLick(1);
        spike_times_psth{n_trial,1} = spk_times_tmp{i_trial} - time_1stLick;
    end
    [psth t_lick] = func_getPSTH_smallBin(spike_times_psth,-3.5,2);
    psth_yes = psth;
    
    subplot(5,6,n_cell); hold on
    plot(t_lick,psth,'b');
    
    
    % no trials, correct
    spike_times_psth = {};
    spk_times_tmp = spk_times_no_correct_all{i_cell,1};
    lick_times_tmp = lick_times_all{i_cell,2};
    n_trial = 0;
    for i_trial = 1:size(spk_times_tmp,1)
        n_trial = n_trial+1;
        time_1stLick = lick_times_tmp{i_trial};
        time_1stLick = time_1stLick(time_1stLick(:,2)>0,2);
        time_1stLick = time_1stLick(1);
        spike_times_psth{n_trial,1} = spk_times_tmp{i_trial} - time_1stLick;
    end
    [psth t_lick] = func_getPSTH_smallBin(spike_times_psth,-3.5,2);
    psth_no = psth;
    
    subplot(5,6,n_cell); hold on
    plot(t_lick,psth,'r');
    xlim([-.2 .8])
    
    ylim_tmp = get(gca,'ylim');
    line([0 0],ylim_tmp,'color','k')
            
end




%% fraction of selectivity neurons
CCF_x = CCF_x_y_slice(:,1)

fract1=[];  fract2=[];  fract3=[];
figure;
for CCF_x_coord = [340:20:480]
    
    if CCF_x_coord == 360
        i_location  = CCF_x<(CCF_x_coord+20) & CCF_x_y_slice(:,3)<320;
    elseif CCF_x_coord == 480
        i_location  = CCF_x>(CCF_x_coord-20) & CCF_x_y_slice(:,3)<320;
    else
        i_location  = CCF_x<(CCF_x_coord+20) & CCF_x>(CCF_x_coord-20) & CCF_x_y_slice(:,3)<320;
    end
    
    fract1(end+1,1) = sum(p_coeff_alpha(i_location,5)<0.01)/sum(~isnan(p_coeff_alpha(i_location,5)));
    fract2(end+1,1) = sum(p_coeff_beta(i_location,5)<0.01)/sum(~isnan(p_coeff_beta(i_location,5)));
    fract3(end+1,1) = sum(Lick_spike_count(i_location,6)<(.01))/sum(~isnan(Lick_spike_count(i_location,6)));
    
end
subplot(3,1,1); hold on; h=bar([340:20:480],fract1); set(h,'facecolor','g'); title('stimulus selective neurons')
subplot(3,1,2); hold on; h=bar([340:20:480],fract2); set(h,'facecolor','m'); title('choice selective neurons')
subplot(3,1,3); hold on; h=bar([340:20:480],fract3); set(h,'facecolor','y'); title('licking movement neurons')
xlabel('CCF ML coordinate')
ylabel('Fraction of selective neurons')



keyboard

% ---------- stat test by bootstrap ----------
peak_stim_btstrp = [];
peak_choice_btstrp = [];
peak_mov_btstrp = [];
%figure
for i_btstrp = 1:1e4
    
    if rem(i_btstrp,1000)==0
        i_btstrp
        %keyboard
    end
    
    i_sample = randsample(size(CCF_x_y_slice,1),size(CCF_x_y_slice,1),'true');
    
    CCF_x_iBtstrp = CCF_x_y_slice(i_sample,1);
    CCF_z_iBtstrp = CCF_x_y_slice(i_sample,3);
    p_coeff_alpha_iBtstrp = p_coeff_alpha(i_sample,:);
    p_coeff_beta_iBtstrp = p_coeff_beta(i_sample,:);
    Lick_spike_count_iBtstrp = Lick_spike_count(i_sample,:);
    
    fract1=[];  fract2=[];  fract3=[];
    for CCF_x_coord = [360:5:480]
        
        i_location  = CCF_x_iBtstrp<(CCF_x_coord+20) & CCF_x_iBtstrp>(CCF_x_coord-20) & CCF_z_iBtstrp<320;
        
        if sum(~isnan(p_coeff_alpha_iBtstrp(i_location,5)))>=50
            fract1(end+1,1) = sum(p_coeff_alpha_iBtstrp(i_location,5)<0.01)/sum(~isnan(p_coeff_alpha_iBtstrp(i_location,5)));
        else
            fract1(end+1,1) = nan;
        end
        if sum(~isnan(p_coeff_beta_iBtstrp(i_location,5)))>=50
            fract2(end+1,1) = sum(p_coeff_beta_iBtstrp(i_location,5)<0.01)/sum(~isnan(p_coeff_beta_iBtstrp(i_location,5)));
        else
            fract2(end+1,1) = nan;
        end
        if sum(~isnan(Lick_spike_count_iBtstrp(i_location,6)))>=50
            fract3(end+1,1) = sum(Lick_spike_count_iBtstrp(i_location,6)<(.01))/sum(~isnan(Lick_spike_count_iBtstrp(i_location,6)));
        else
            fract3(end+1,1) = nan;
        end        
    end
    
    CCF_x_coord = [360:5:480];
    
    [dummy i_peak1] = max(fract1);
    peak_stim_btstrp(i_btstrp,1) = CCF_x_coord(i_peak1);
    
    [dummy i_peak2] = max(fract2);
    peak_choice_btstrp(i_btstrp,1) = CCF_x_coord(i_peak2);
    
    [dummy i_peak3] = max(fract3);
    peak_mov_btstrp(i_btstrp,1) = CCF_x_coord(i_peak3);
    
    %hold on; plot([360:5:480],fract1,'g'); plot(CCF_x_coord(i_peak1),fract1(i_peak1),'ok');
    %hold on; plot([360:5:480],fract2,'m'); plot(CCF_x_coord(i_peak2),fract2(i_peak2),'ok');
    %hold on; plot([360:5:480],fract3,'y'); plot(CCF_x_coord(i_peak3),fract3(i_peak3),'ok');

end
disp('============ Stim vs. Lick Movement ========')
p = sum((peak_stim_btstrp-peak_mov_btstrp)<0)/size(peak_stim_btstrp,1);
disp(['p=',num2str(p)]);

disp('============ Choice vs. Lick Movement ========')
p = sum((peak_choice_btstrp-peak_mov_btstrp)<0)/size(peak_stim_btstrp,1);
disp(['p=',num2str(p)]);

disp('============ Stim vs. Choice ========')
p = sum((peak_stim_btstrp-peak_choice_btstrp)<0)/size(peak_stim_btstrp,1);
disp(['p=',num2str(p)]);
