%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 2f SC selectivity map


clear all
close all

addpath('../func/')

load Figure_2_f_data

sig_selective(isnan(sig_selective))=1;
sig_selective = sig_selective<.01; % <----------- This p value can be adjusted


%% ---------------- view contra/ipsi selectivity in CCF space --------------------------

% sample epoch
% load CCF images
CCF_sections_all = 860:20:920;

% figure;
for i_section = 1:length(CCF_sections_all)    
    im_CCF = imread(['..\CCF_alignment\ARA_template\AllenRefVolCoronal_100',num2str(CCF_sections_all(i_section)),'.tif']);
end

%% stimulus, choice, movement selectivity, by epoch
p_coeff_alpha = []; 
p_coeff_beta = []; 
p_coeff_sigma = []; 
for i_unit = 1:size(sig_selective,1)
    
    spike_times1 = spike_times_all{i_unit,1};       %[yes no yes_error no_error]
    spike_times2 = spike_times_all{i_unit,2};       %[yes no yes_error no_error]
    spike_times3 = spike_times_all{i_unit,3};       %[yes no yes_error no_error]
    spike_times4= spike_times_all{i_unit,4};       %[yes no yes_error no_error]

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

%% Summary Stimulus, Choice, Movement selectivity MAP (for these maps, all neurons from SC are collapsed together, see detailed Maps at different sections).
% load CCF images
im_CCF = imread(['..\CCF_alignment\ARA_template\AllenRefVolCoronal_100880.tif']);

im_CCF_anatomy = imread(['.\ALM2x_mergedRGB_slice880_cltr1247_P55_NL09.tif']);
im_CCF_anatomy = imresize(im_CCF_anatomy,2);

figure;
for i_tmp = 1:3
    subplot(3,1,i_tmp);
    image(im_CCF); colormap gray(250); hold on
    plot(Unit_CCF(:,2),Unit_CCF(:,3),'ow','markersize',1)
    xlim([330 520]); ylim([190 340]);
end

% stimulus selectivity
%max_sel = max(max(abs(-log10([p_coeff_alpha(:,5)]))));
selectivity_stim = abs(spike_count_yes_stim - spike_count_no_stim);
my_dot_size = (selectivity_stim(:,5))/10*8+2;     % selectity amplitude
my_dot_size(my_dot_size>10)=10;             % cap at 10
for i_unit = 1:size(Unit_CCF,1)
    
    % sample + delay
    if p_coeff_alpha(i_unit,5)<0.01
        subplot(3,1,1);
        %plot(Unit_CCF(i_unit,2),Unit_CCF(i_unit,3),'og','markerfacecolor','g','markersize', -log10(p_coeff_alpha(i_unit,5))/max_sel*20);
        plot(Unit_CCF(i_unit,2),Unit_CCF(i_unit,3),'og','markerfacecolor','g','markersize', my_dot_size(i_unit,1));
    end
    
end
subplot(3,1,1); title('Stimulus selectivity, Sample + Delay')

% choice selectivity
selectivity_choice = abs(spike_count_lickYes_choice - spike_count_lickNo_choice);
my_dot_size = (selectivity_choice(:,5))/10*8+2;     % selectity amplitude
my_dot_size(my_dot_size>10)=10;             % cap at 10
for i_unit = 1:size(Unit_CCF,1)
    
    % sample + delay
    if p_coeff_beta(i_unit,5)<0.01
        subplot(3,1,2);
        plot(Unit_CCF(i_unit,2),Unit_CCF(i_unit,3),'om','markerfacecolor','m','markersize', my_dot_size(i_unit,1));
    end
    
end
subplot(3,1,2); title('Choice selectivity, Sample + Delay')


% Licking modulation
Mod_all = abs(Lick_spike_count(:,5));

my_dot_size = (Mod_all)*10;%/1*5+2;     % selectity amplitude
my_dot_size(my_dot_size>10)=10;             % cap at 10
for i_unit = 1:size(Lick_spike_count,1)
    
    if Lick_spike_count(i_unit,6)<(.01/size(Lick_spike_count,1))
        subplot(3,1,3);
        plot(Unit_CCF(i_unit,2),Unit_CCF(i_unit,3),'oy','markerfacecolor','y','markersize', my_dot_size(i_unit,1));
    end
    
end
title('Licking modulation, Response')




%% fraction of selective neurons
CCF_x = Unit_CCF(:,2);

fract1=[];  fract2=[];  fract3=[];
figure;
for CCF_x_coord = [380:20:480]
    
    if CCF_x_coord == 360
        i_location  = CCF_x<(CCF_x_coord+20);
    elseif CCF_x_coord == 480
        i_location  = CCF_x>(CCF_x_coord-20);
    else
        i_location  = CCF_x<(CCF_x_coord+20) & CCF_x>(CCF_x_coord-20);
    end
    
    fract1(end+1,1) = sum(p_coeff_alpha(i_location,5)<0.01)/sum(~isnan(p_coeff_alpha(i_location,5)));
    fract2(end+1,1) = sum(p_coeff_beta(i_location,5)<0.01)/sum(~isnan(p_coeff_beta(i_location,5)));
    fract3(end+1,1) = sum(Lick_spike_count(i_location,6)<(.01))/sum(~isnan(Lick_spike_count(i_location,6)));
    
end
subplot(3,1,1); hold on; h=bar([380:20:480],fract1); set(h,'facecolor','g'); title('stimulus selective neurons')
subplot(3,1,2); hold on; h=bar([380:20:480],fract2); set(h,'facecolor','m'); title('stimulus selective neurons')
subplot(3,1,3); hold on; h=bar([380:20:480],fract3); set(h,'facecolor','y'); title('stimulus selective neurons')
xlabel('CCF ML coordinate')
ylabel('Fraction of selective neurons')
