%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 7d direct photoinhibition of SC GABAergic neurons

% flexed ArchT in left SC of gad2- or vgat-cre mice, left SC fiber implantation

clear all
close all

load Figure_7_d_data

%% check if the data seems right
% hist(StimOnTime_allSession(Sample_Delay_allSession==1 & LickEarly_allSession==0 & AOM_data_allSession>0))
% 
% unique(StimDur_allSession(LickEarly_allSession==0))

%%
perf_Raw_all = [];

n_mice = 0;
n_plot = 0;
for i_mice = unique(Session_type_allSession(:,1))'
    
    n_mice = n_mice+1;
    n_plot = n_plot+1;
    
    if n_plot>8
        figure
        n_plot = 1;
    end
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
    

    figure(1);
    subplot(2,4,n_plot); hold on
    errorbar(1:4,Y_perf(:,1),Y_perf(:,2),'color','b')
    errorbar(1:4,Y_perf(:,3),Y_perf(:,4),'color','r')
    for i_tmp = 1:4
        text(i_tmp,Y_perf(i_tmp,1)+.1,num2str(n_trials(i_tmp,1)));
        text(i_tmp,Y_perf(i_tmp,3)+.1,num2str(n_trials(i_tmp,2)));
    end
    ylabel('Fraction correct')
    xlim([0 5])
    set(gca, 'XTick', [1 2 3 4])
    set(gca,'XTickLabel',str2mat('1','2','3','4'))
    set(gca, 'XTickLabel', {['Ctrl'] ['S'] ['D'] ['R'] })
    ylim([0 1])
    title(['Mouse ',num2str(i_mice)]);
    
    
    % save data
    perf_Raw_all(:,1,n_mice) = Y_perf(:,1);
    perf_Raw_all(:,2,n_mice) = Y_perf(:,3);
        
    
    
end



% plot performance for all mice
    perf_yes = squeeze(perf_Raw_all(:,1,:));
    perf_no = squeeze(perf_Raw_all(:,2,:));

    figure; hold on
    plot(1:4,perf_no,'color',[1 .6 .6]);
    plot(1:4,mean(perf_no,2),'-or','linewidth',2,'markerfacecolor','r');

    plot(1:4,perf_yes,'color',[.6 .6 1]);
    plot(1:4,mean(perf_yes,2),'-ob','linewidth',2,'markerfacecolor','b');
    title('All mice')
    ylim([0.4 1])
    xlim([0 5])
    set(gca, 'XTick', [1 2 3 4])
    set(gca,'XTickLabel',str2mat('1','2','3','4'))
    set(gca, 'XTickLabel', {['Ctrl'] ['S'] ['D'] ['R'] })
    pbaspect([1 1.5 1]) 
    set(findall(gcf,'-property','FontSize'),'FontSize',16)






%% ============ statistics using bootsrap ===============

i_sel_mice = unique(Session_type_allSession(:,1))'


delta_perf_btstrp_contra = [];
delta_perf_btstrp_ipsi = [];

for i_btstrp = 1:10^4
    if rem(i_btstrp,1000)==0
        disp(i_btstrp)
    end
    
    i_sample_animal = randsample(i_sel_mice,length(i_sel_mice),'true');
    
    Y_perf_control = [];
    Y_perf_stim = [];
    
    n_mice = 0;
    for i_mice = i_sample_animal
               
        n_mice = n_mice+1;
        
        AOM_data_tmp = [];
        StimTrials_tmp = [];
        Session_type_tmp = [];
        LickEarly_tmp = [];
        Sample_Delay_tmp = [];
        StimDur_tmp = [];
        R_hit_tmp = [];
        R_miss_tmp = [];
        L_hit_tmp = [];
        L_miss_tmp = [];
        StimOnTime_tmp = [];
        
        i_selectSession = unique(Session_type_allSession(Session_type_allSession(:,1)==i_mice, 2));
        n_session = size(i_selectSession,1);
        i_sample_session = randsample(n_session,n_session,'true');
        i_sample_session = i_selectSession(i_sample_session);
        
        
        for i_session = i_sample_session'
            
            i_selectTrial = find(Session_type_allSession(:,1)==i_mice & Session_type_allSession(:,2)==i_session);
            n_trials = size(i_selectTrial,1);
            i_sample_trialID = randsample(n_trials,n_trials,'true');
            i_sample_trialID = i_selectTrial(i_sample_trialID,1);
            
            
            AOM_data_tmp        = cat(1, AOM_data_tmp, AOM_data_allSession(i_sample_trialID,:));
            StimTrials_tmp      = cat(1, StimTrials_tmp, StimTrials_allSession(i_sample_trialID,:));
            Session_type_tmp    = cat(1, Session_type_tmp, Session_type_allSession(i_sample_trialID,:));
            LickEarly_tmp       = cat(1, LickEarly_tmp, LickEarly_allSession(i_sample_trialID,:));
            Sample_Delay_tmp    = cat(1, Sample_Delay_tmp, Sample_Delay_allSession(i_sample_trialID,:));
            StimDur_tmp         = cat(1, StimDur_tmp, StimDur_allSession(i_sample_trialID,:));
            R_hit_tmp           = cat(1, R_hit_tmp, R_hit_allSession(i_sample_trialID,:));
            R_miss_tmp          = cat(1, R_miss_tmp, R_miss_allSession(i_sample_trialID,:));
            L_hit_tmp           = cat(1, L_hit_tmp, L_hit_allSession(i_sample_trialID,:));
            L_miss_tmp          = cat(1, L_miss_tmp, L_miss_allSession(i_sample_trialID,:));
            StimOnTime_tmp      = cat(1, StimOnTime_tmp, StimOnTime_allSession(i_sample_trialID,:));
            
        end
            
        % comput performance
        % invert the trial type based on stimulated SC
        i_select_leftSC = find(Session_type_tmp(:,3)==1);
%         i_select_rightSC = find(Session_type_tmp(:,4)==2);
        
        Contra_hit_tmp(i_select_leftSC,:)  = R_hit_tmp(i_select_leftSC);
%         Contra_hit_tmp(i_select_rightSC,:) = L_hit_tmp(i_select_rightSC);
        Contra_miss_tmp(i_select_leftSC,:)  = R_miss_tmp(i_select_leftSC);
%         Contra_miss_tmp(i_select_rightSC,:) = L_miss_tmp(i_select_rightSC);
        
        Ipsi_hit_tmp(i_select_leftSC,:)  = L_hit_tmp(i_select_leftSC);
%         Ipsi_hit_tmp(i_select_rightSC,:) = R_hit_tmp(i_select_rightSC);
        Ipsi_miss_tmp(i_select_leftSC,:)  = L_miss_tmp(i_select_leftSC);
%         Ipsi_miss_tmp(i_select_rightSC,:) = R_miss_tmp(i_select_rightSC);
        
        
        % control
        X_type = [];
        Y_perf = [];
        n_trials = [];
        
        X_type(end+1,1) = size(X_type,1)+1;
        
        i_select = find(AOM_data_tmp == 0 & StimTrials_tmp >0 & LickEarly_tmp == 0 & Session_type_tmp(:,1)==i_mice);
        
        % Contra Trials
        Y_perf(end+1,1) = sum(Contra_hit_tmp(i_select))/sum(Contra_hit_tmp(i_select)|Contra_miss_tmp(i_select));
        
        % Ipsi Trials
        Y_perf(end,2) = sum(Ipsi_hit_tmp(i_select))/sum(Ipsi_hit_tmp(i_select)|Ipsi_miss_tmp(i_select));
        
        
        % stim trials
        for i_condition = 1
            
            X_type(end+1,1) = size(X_type,1)+1;
            
            if i_condition==1
             
                i_select = (Sample_Delay_tmp==2 & AOM_data_tmp>1 & AOM_data_tmp<6 & LickEarly_tmp==0 & StimTrials_tmp >0 & Session_type_tmp(:,1)==i_mice & StimDur_tmp==1.3);
                 end
            
            % Contra Trials
            Y_perf(end+1,1) = sum(Contra_hit_tmp(i_select))/sum(Contra_hit_tmp(i_select)|Contra_miss_tmp(i_select));
            
            % Ipsi Trials
            Y_perf(end,2) = sum(Ipsi_hit_tmp(i_select))/sum(Ipsi_hit_tmp(i_select)|Ipsi_miss_tmp(i_select));

        end
        
        perf_Raw_btstrp(:,1,n_mice) = Y_perf(:,1);
        perf_Raw_btstrp(:,2,n_mice) = Y_perf(:,2);

    end
    
    perf_tmp = mean(perf_Raw_btstrp,3);

    delta_perf_btstrp_contra(i_btstrp,:) = perf_tmp(:,1)-perf_tmp(1,1);
    delta_perf_btstrp_ipsi(i_btstrp,:) = perf_tmp(:,2)-perf_tmp(1,2);
  
end
delta_perf_btstrp_contra(isnan(sum(delta_perf_btstrp_contra,2)),:) = [];
delta_perf_btstrp_ipsi(isnan(sum(delta_perf_btstrp_ipsi,2)),:) = [];

disp('============= delay perf Contra, p value ==============')
p=sum(delta_perf_btstrp_contra(:,2:end)<=0)/size(delta_perf_btstrp_contra,1)

disp('=============delay perf ipsi, p value ==============')
p=sum(delta_perf_btstrp_ipsi(:,2:end)>=0)/size(delta_perf_btstrp_ipsi,1)


