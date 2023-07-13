%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 2e ALM selectivity map

clear all
close all
addpath('../func/');

load Figure_2_e_data


im_CCF = imread(['..\CCF_alignment\ARA_template\AllenRefVolCoronal_100275.tif']);

im_CCF_anatomy = imread(['.\ALM2x_mergedRGB_slice300_cltr1247_P55_NL09.tif']);
im_CCF_anatomy = imresize(im_CCF_anatomy,2);

%% Summary Stimulus, Choice, Movement selectivity MAP (for these maps, the most posterior neurons are exclude because they are outside of ALM and including them make them appear to be outside the brain).
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

% Licking Movement activity
lick_sel_cell = Lick_spike_count(:,6)<(.01/size(Lick_spike_count,1));

Mod_all = abs(Lick_spike_count(:,5));


% Licking modulation
for i_unit = 1:size(Lick_spike_count,1)
    
    if Lick_spike_count(i_unit,6)<(.01/size(Lick_spike_count,1)) & CCF_x_y_slice(i_unit,3)<320
        subplot(3,1,3);
        plot(CCF_x_y_slice(i_unit,1),CCF_x_y_slice(i_unit,2),'oy','markerfacecolor','y','markersize', Mod_all(i_unit)*20);
    end
    
end
title('Licking modulation, Response')



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

