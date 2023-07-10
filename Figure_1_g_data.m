%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 1g, ignore rate

clear all
close all

load Figure_1_g_IgnoreRate



%% individual mouse
perf_Raw_all = [];

n_mice = 0;
for i_mice = unique(Session_type_allSession(:,1))'
    
    
    n_mice = n_mice+1;
    
    %% ------------ plot the performance --------------
    
    
    % control
    X_type = [];
    Y_perf = [];
    n_trials = [];
    
    X_type(end+1,1) = size(X_type,1)+1;
    
    i_select = find(AOM_data_allSession == 0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice);
    
    
    % Ignore rate, alls Trials
    Y_perf(end+1,1) = sum(R_ignore_allSession(i_select)|L_ignore_allSession(i_select))/length(i_select);
    for i=1:1000
        tmp = R_ignore_allSession(i_select);
        perf_tmp(i) = mean(tmp(randsample(length(tmp),length(tmp),1)));
    end
    Y_perf(end,2) = std(perf_tmp);
    n_trials(end+1,1) = length(i_select);

    
    
    % stim trials
    for i_condition = 1:3
        
        X_type(end+1,1) = size(X_type,1)+1;
        
        if i_condition==1
            i_select = find(Sample_Delay_allSession==1 & AOM_data_allSession>1.5 & LickEarly_allSession==0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice & StimDur_allSession>0.6);
        elseif i_condition==2
            i_select = find(Sample_Delay_allSession==2 & AOM_data_allSession>1.5 & LickEarly_allSession==0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice & StimDur_allSession>0.6);
        elseif i_condition==3
            i_select = find(Sample_Delay_allSession==3 & AOM_data_allSession>1.5 & LickEarly_allSession==0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice & StimDur_allSession>0.6);
        end
        
        % Ignore rate, alls Trials
        Y_perf(end+1,1) = sum(R_ignore_allSession(i_select)|L_ignore_allSession(i_select))/length(i_select);
        for i=1:1000
            tmp = R_ignore_allSession(i_select);
            perf_tmp(i) = mean(tmp(randsample(length(tmp),length(tmp),1)));
        end
        Y_perf(end,2) = std(perf_tmp);
        n_trials(end+1,1) = length(i_select);
        
    end
    

    figure(1);
    subplot(2,3,i_mice); hold on
    errorbar(1:4,Y_perf(:,1),Y_perf(:,2),'color','k')
    for i_tmp = 1:4
        text(i_tmp,Y_perf(i_tmp,1)+.1,num2str(n_trials(i_tmp,1)));
    end
    ylabel('Fraction ignored')
    ylim([0 1])
    title(['Mouse ',num2str(i_mice)]);
    
    
    % save data
    perf_Raw_all(:,1,n_mice) = Y_perf(:,1);
    
    
end




% plot performance for all mice
perf_all = squeeze(perf_Raw_all(:,1,:));

figure; hold on
plot(1:4,perf_all,'color',[.6 .6 .6]);
plot(1:4,mean(perf_all,2),'-ok','linewidth',2,'markerfacecolor','k');
line([0 5],[.5 .5],'color','k','lineStyle',':')
ylim([0 1])
xlim([0 5])
ylabel('Fraction ignored')
xlabel('Control  Sample  Delay  Response')

