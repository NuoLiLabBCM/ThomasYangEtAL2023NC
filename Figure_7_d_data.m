%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 7d direct photoinhibition of SC GABAergic neurons

% flex-ArchT in left SC of gad2- or vgat-cre mice, left SC fiber implantation

clear all
close all

load Figure_7_d_data

perf_Raw_all = [];

n_mice = 0;
n_plot = 0;
for i_mice = unique(Session_type_allSession(:,1))'
    
    n_mice = n_mice+1;
    n_plot = n_plot+1;

    %% ------------ plot the performance --------------
    
    
    % control
    X_type = [];
    Y_perf = [];
    n_trials = [];
    
    X_type(end+1,1) = size(X_type,1)+1;
    
    i_select = find(AOM_data_allSession == 0 & StimTrials_allSession >0 & LickEarly_allSession == 0 & Session_type_allSession(:,1)==i_mice);
    
    
    % Yes Trials
    Y_perf(end+1,1) = sum(R_hit_allSession(i_select))/sum(R_hit_allSession(i_select)|R_miss_allSession(i_select));
    for i=1:1000
        tmp = R_hit_allSession(i_select);
        i_tmp = find(R_hit_allSession(i_select)|R_miss_allSession(i_select));
        tmp = tmp(i_tmp);
        perf_tmp(i) = mean(tmp(randsample(length(tmp),length(tmp),1)));
    end
    Y_perf(end,2) = std(perf_tmp);
    n_trials(end+1,1) = sum(R_hit_allSession(i_select)|R_miss_allSession(i_select));
    
    % No Trials
    Y_perf(end,3) = sum(L_hit_allSession(i_select))/sum(L_hit_allSession(i_select)|L_miss_allSession(i_select));
    for i=1:1000
        tmp = L_hit_allSession(i_select);
        i_tmp = find(L_hit_allSession(i_select)|L_miss_allSession(i_select));
        tmp = tmp(i_tmp);
        perf_tmp(i) = mean(tmp(randsample(length(tmp),length(tmp),1)));
    end
    Y_perf(end,4) = std(perf_tmp);
    n_trials(end,2) = sum(L_hit_allSession(i_select)|L_miss_allSession(i_select));    
    
    
    % stim trials
    for i_condition = 1:3
        
        X_type(end+1,1) = size(X_type,1)+1;
        
        if i_condition==1
            i_select = (Sample_Delay_allSession==1 & AOM_data_allSession>0 & AOM_data_allSession<6 & LickEarly_allSession==0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice & StimDur_allSession==1.3);
        elseif i_condition==2
            i_select = (Sample_Delay_allSession==2 & AOM_data_allSession>0 & AOM_data_allSession<6 & LickEarly_allSession==0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice & StimDur_allSession==1.3);
        elseif i_condition==3
            i_select = (Sample_Delay_allSession==3 & AOM_data_allSession>0 & AOM_data_allSession<6 & LickEarly_allSession==0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice & StimDur_allSession==1.3);
        end
        
        % Yes Trials
        Y_perf(end+1,1) = sum(R_hit_allSession(i_select))/sum(R_hit_allSession(i_select)|R_miss_allSession(i_select));
        for i=1:1000
            tmp = R_hit_allSession(i_select);
            i_tmp = find(R_hit_allSession(i_select)|R_miss_allSession(i_select));
            tmp = tmp(i_tmp);
            perf_tmp(i) = mean(tmp(randsample(length(tmp),length(tmp),1)));
        end
        Y_perf(end,2) = std(perf_tmp);
        n_trials(end+1,1) = sum(R_hit_allSession(i_select)|R_miss_allSession(i_select));
        
        % No Trials
        Y_perf(end,3) = sum(L_hit_allSession(i_select))/sum(L_hit_allSession(i_select)|L_miss_allSession(i_select));
        for i=1:1000
            tmp = L_hit_allSession(i_select);
            i_tmp = find(L_hit_allSession(i_select)|L_miss_allSession(i_select));
            tmp = tmp(i_tmp);
            perf_tmp(i) = mean(tmp(randsample(length(tmp),length(tmp),1)));
        end
        Y_perf(end,4) = std(perf_tmp);
        n_trials(end,2) = sum(L_hit_allSession(i_select)|L_miss_allSession(i_select));
        
        
    end
    
     % save data
    perf_Raw_all(:,1,n_mice) = Y_perf(:,1);
    perf_Raw_all(:,2,n_mice) = Y_perf(:,3);
        
    
    
end


% plot performance for all mice
    perf_yes = squeeze(perf_Raw_all(:,1,:));
    perf_no = squeeze(perf_Raw_all(:,2,:));

    figure; hold on
    plot(1:2,perf_no([1 3],:),'color',[1 .6 .6]);
    plot(1:2,mean(perf_no([1 3],:),2),'-or','linewidth',2,'markerfacecolor','r');

    plot(1:2,perf_yes([1 3],:),'color',[.6 .6 1]);
    plot(1:2,mean(perf_yes([1 3],:),2),'-ob','linewidth',2,'markerfacecolor','b');
    title('All mice')
    ylim([0.4 1])
    xlim([0.5 2.5])
    ylabel('Fraction correct');
    set(gca, 'XTick', [1 2])
    set(gca,'XTickLabel',str2mat('1','2'))
    set(gca, 'XTickLabel', {['No stim'] ['D']})
    set(findall(gcf,'-property','FontSize'),'FontSize',16)    
