%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 3a ALM selectivity map

clear all
close all
addpath('../func/');

load Figure_3_a_data


im_CCF = imread(['..\CCF_alignment\ARA_template\AllenRefVolCoronal_100275.tif']);

im_CCF_anatomy = imread(['.\ALM2x_mergedRGB_slice300_cltr1247_P55_NL09.tif']);
im_CCF_anatomy = imresize(im_CCF_anatomy,2);


% plot contra vs. ipsi selectivity 
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







%% Fig. 3a SC selectivity map
% SC selectivity map, contra vs. ipsi plot, collapsing all cells onto one section
% load CCF images
clear all

im_CCF = imread(['..\CCF_alignment\ARA_template\AllenRefVolCoronal_100880.tif']);

im_CCF_anatomy = imread(['.\ALM2x_mergedRGB_slice880_cltr1247_P55_NL09.tif']);
im_CCF_anatomy = imresize(im_CCF_anatomy,2);

load Figure_2_f_data

sig_selective(isnan(sig_selective))=1;
sig_selective = sig_selective<.01; % <----------- This p value can be adjusted
selectivity_all = spk_count_yes_all-spk_count_no_all;   %[sample delay response whole_trial baseline]

figure
for i_epoch = 1:3
    subplot(1,3,i_epoch);
    imagesc(im_CCF_anatomy);
    colormap('gray')
    hold on    
    plot(Unit_CCF(:,2),Unit_CCF(:,3),'ow','markerfacecolor','w','markersize',1.5)
    xlim([330 520]); ylim([190 340]);

    
    for i_unit = 1:size(sig_selective,1)
        
        % set color and dot size based on selectivity
        if selectivity_all(i_unit,i_epoch)>0
            my_color = 'b';     % contra
        else
            my_color = 'r';     % ipsi
        end
        
        my_dot_size = (abs(selectivity_all(i_unit,i_epoch))/20)*6+2;     % selectity amplitude
        my_dot_size(my_dot_size>10)=10;             % cap at 10
        
        % plot unit in CCF
        if sig_selective(i_unit,i_epoch)
            plot(Unit_CCF(i_unit,2),Unit_CCF(i_unit,3),'ow','markerfacecolor',my_color,'markersize',my_dot_size);
        end
    end
end


