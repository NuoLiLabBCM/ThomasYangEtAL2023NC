%% Thomas_Yang_et al, 2023 @Nuo Li lab
%% Fig. 1h, performance with unilateral SC photoinhibition with vgat mice

clear all
close all

load Figure_1_h_data


for i_mice = unique(Session_type_allSession(:,1))'
    disp(['======= mouse ',num2str(i_mice),' [leftSC rightSC] ==========']);
    [sum(Session_type_allSession(Session_type_allSession(:,1)==i_mice,4)==1) sum(Session_type_allSession(Session_type_allSession(:,1)==i_mice,4)==2)]
    disp(['power']);
    unique(AOM_data_allSession(Session_type_allSession(:,1)==i_mice))
end

SC_side =[  % 1--left;  2--right  for the last 2 mice, only data from 1 hemisphere will be used
    2;...
    2;...
    2;...
    1;...
    1;...
    1;...
    2;...
    1;...
    ];

power_thresh = [   % power level used, for each mouse, only use the strongest powers level tested
    1;...
    1;...
    1;...
    .5;...
    .6;...
    1;...
    .6;...
    .6;...
    ];

%% performance, control vs 4 epochs
perf_Raw_all = [];

n_mice = 0;
for i_mice = unique(Session_type_allSession(:,1))'
    
    n_mice = n_mice+1;
    
    
    %% ------------ plot the performance --------------
    
    % invert the trial type based on stimulated SC
    i_select_leftSC = find(Session_type_allSession(:,4)==1);
    i_select_rightSC = find(Session_type_allSession(:,4)==2);
    
    Contra_hit_allSession(i_select_leftSC,:)  = R_hit_allSession(i_select_leftSC);
    Contra_hit_allSession(i_select_rightSC,:) = L_hit_allSession(i_select_rightSC);
    Contra_miss_allSession(i_select_leftSC,:)  = R_miss_allSession(i_select_leftSC);
    Contra_miss_allSession(i_select_rightSC,:) = L_miss_allSession(i_select_rightSC);
    Contra_ignore_allSession(i_select_leftSC,:)  = R_ignore_allSession(i_select_leftSC);
    Contra_ignore_allSession(i_select_rightSC,:) = L_ignore_allSession(i_select_rightSC);
    
    Ipsi_hit_allSession(i_select_leftSC,:)  = L_hit_allSession(i_select_leftSC);
    Ipsi_hit_allSession(i_select_rightSC,:) = R_hit_allSession(i_select_rightSC);
    Ipsi_miss_allSession(i_select_leftSC,:)  = L_miss_allSession(i_select_leftSC);
    Ipsi_miss_allSession(i_select_rightSC,:) = R_miss_allSession(i_select_rightSC);
    Ipsi_ignore_allSession(i_select_leftSC,:)  = L_ignore_allSession(i_select_leftSC);
    Ipsi_ignore_allSession(i_select_rightSC,:) = R_ignore_allSession(i_select_rightSC);
    
    
    % control
    X_type = [];
    Y_perf = [];
    n_trials = [];
    
    X_type(end+1,1) = size(X_type,1)+1;
    
    i_select = find(AOM_data_allSession == 0 & StimTrials_allSession >0 & LickEarly_allSession == 0 & Session_type_allSession(:,1)==i_mice & Session_type_allSession(:,4)==SC_side(n_mice));
    
    
    % Contra Trials
    Y_perf(end+1,1) = sum(Contra_hit_allSession(i_select))/sum(Contra_hit_allSession(i_select)|Contra_miss_allSession(i_select));
    for i=1:1000
        tmp = Contra_hit_allSession(i_select);
        i_tmp = find(Contra_hit_allSession(i_select)|Contra_miss_allSession(i_select));
        tmp = tmp(i_tmp);
        perf_tmp(i) = mean(tmp(randsample(length(tmp),length(tmp),1)));
    end
    Y_perf(end,2) = std(perf_tmp);
    n_trials(end+1,1) = sum(Contra_hit_allSession(i_select)|Contra_miss_allSession(i_select));
    
    % Ipsi Trials
    Y_perf(end,3) = sum(Ipsi_hit_allSession(i_select))/sum(Ipsi_hit_allSession(i_select)|Ipsi_miss_allSession(i_select));
    for i=1:1000
        tmp = Ipsi_hit_allSession(i_select);
        i_tmp = find(Ipsi_hit_allSession(i_select)|Ipsi_miss_allSession(i_select));
        tmp = tmp(i_tmp);
        perf_tmp(i) = mean(tmp(randsample(length(tmp),length(tmp),1)));
    end
    Y_perf(end,4) = std(perf_tmp);
    n_trials(end,2) = sum(Ipsi_hit_allSession(i_select)|Ipsi_miss_allSession(i_select));
    
    
    
    % stim trials
    for i_condition = 1:4
        
        X_type(end+1,1) = size(X_type,1)+1;
        
        if i_condition==1
            i_select = (Sample_Delay_allSession==1 & AOM_data_allSession>=power_thresh(n_mice)  & StimOnTime_allSession < 1 & LickEarly_allSession==0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice & StimDur_allSession>0 & Session_type_allSession(:,4)==SC_side(n_mice));
        elseif i_condition==2
            i_select = (Sample_Delay_allSession==1 & AOM_data_allSession>=power_thresh(n_mice)  & StimOnTime_allSession > 1 & LickEarly_allSession==0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice & StimDur_allSession>0 & Session_type_allSession(:,4)==SC_side(n_mice));
        elseif i_condition==3
            i_select = (Sample_Delay_allSession==2 & AOM_data_allSession>=power_thresh(n_mice)  & StimOnTime_allSession < 2 & LickEarly_allSession==0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice & StimDur_allSession>0 & Session_type_allSession(:,4)==SC_side(n_mice));
        elseif i_condition==4
            i_select = (Sample_Delay_allSession==2 & AOM_data_allSession>=power_thresh(n_mice)  & StimOnTime_allSession > 2 & LickEarly_allSession==0 & StimTrials_allSession >0 & Session_type_allSession(:,1)==i_mice & StimDur_allSession>0 & Session_type_allSession(:,4)==SC_side(n_mice));
        end
        
        % Contra Trials
        Y_perf(end+1,1) = sum(Contra_hit_allSession(i_select))/sum(Contra_hit_allSession(i_select)|Contra_miss_allSession(i_select));
        for i=1:1000
            tmp = Contra_hit_allSession(i_select);
            i_tmp = find(Contra_hit_allSession(i_select)|Contra_miss_allSession(i_select));
            tmp = tmp(i_tmp);
            perf_tmp(i) = mean(tmp(randsample(length(tmp),length(tmp),1)));
        end
        Y_perf(end,2) = std(perf_tmp);
        n_trials(end+1,1) = sum(Contra_hit_allSession(i_select)|Contra_miss_allSession(i_select));
        
        % Ipsi Trials
        Y_perf(end,3) = sum(Ipsi_hit_allSession(i_select))/sum(Ipsi_hit_allSession(i_select)|Ipsi_miss_allSession(i_select));
        for i=1:1000
            tmp = Ipsi_hit_allSession(i_select);
            i_tmp = find(Ipsi_hit_allSession(i_select)|Ipsi_miss_allSession(i_select));
            tmp = tmp(i_tmp);
            perf_tmp(i) = mean(tmp(randsample(length(tmp),length(tmp),1)));
        end
        Y_perf(end,4) = std(perf_tmp);
        n_trials(end,2) = sum(Ipsi_hit_allSession(i_select)|Ipsi_miss_allSession(i_select));
        
    end
    

    % save data
    perf_Raw_all(:,1,n_mice) = Y_perf(:,1);
    perf_Raw_all(:,2,n_mice) = Y_perf(:,3);
    
    
end



% plot performance for all mice
perf_yes = squeeze(perf_Raw_all(:,1,:));
perf_no = squeeze(perf_Raw_all(:,2,:));

figure; hold on
plot(1:5,perf_no,'color',[1 .6 .6]);
plot(1:5,mean(perf_no,2),'-or','linewidth',2,'markerfacecolor','r');

plot(1:5,perf_yes,'color',[.6 .6 1]);
plot(1:5,mean(perf_yes,2),'-ob','linewidth',2,'markerfacecolor','b');
title('All mice')
ylim([0 1])
xlim([0 6])






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
        
        i_selectSession = unique(Session_type_allSession(Session_type_allSession(:,1)==i_mice & Session_type_allSession(:,4)==SC_side(i_mice), 2));
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
        i_select_leftSC = find(Session_type_tmp(:,4)==1);
        i_select_rightSC = find(Session_type_tmp(:,4)==2);
        
        Contra_hit_tmp(i_select_leftSC,:)  = R_hit_tmp(i_select_leftSC);
        Contra_hit_tmp(i_select_rightSC,:) = L_hit_tmp(i_select_rightSC);
        Contra_miss_tmp(i_select_leftSC,:)  = R_miss_tmp(i_select_leftSC);
        Contra_miss_tmp(i_select_rightSC,:) = L_miss_tmp(i_select_rightSC);
        
        Ipsi_hit_tmp(i_select_leftSC,:)  = L_hit_tmp(i_select_leftSC);
        Ipsi_hit_tmp(i_select_rightSC,:) = R_hit_tmp(i_select_rightSC);
        Ipsi_miss_tmp(i_select_leftSC,:)  = L_miss_tmp(i_select_leftSC);
        Ipsi_miss_tmp(i_select_rightSC,:) = R_miss_tmp(i_select_rightSC);
        
        
        % control
        X_type = [];
        Y_perf = [];
        n_trials = [];
        
        X_type(end+1,1) = size(X_type,1)+1;
        
        i_select = find(AOM_data_tmp == 0 & StimTrials_tmp >0 & LickEarly_tmp == 0 & Session_type_tmp(:,1)==i_mice & Session_type_tmp(:,4)==SC_side(i_mice));
        
        % Contra Trials
        Y_perf(end+1,1) = sum(Contra_hit_tmp(i_select))/sum(Contra_hit_tmp(i_select)|Contra_miss_tmp(i_select));
        
        % Ipsi Trials
        Y_perf(end,2) = sum(Ipsi_hit_tmp(i_select))/sum(Ipsi_hit_tmp(i_select)|Ipsi_miss_tmp(i_select));

        
        % stim trials
        for i_condition = 1:4
            
            X_type(end+1,1) = size(X_type,1)+1;
            
            if i_condition==1
                i_select = (Sample_Delay_tmp==1 & AOM_data_tmp>=power_thresh(i_mice)  & StimOnTime_tmp < 1 & LickEarly_tmp==0 & StimTrials_tmp >0 & Session_type_tmp(:,1)==i_mice & StimDur_tmp>0 & Session_type_tmp(:,4)==SC_side(i_mice));
            elseif i_condition==2
                i_select = (Sample_Delay_tmp==1 & AOM_data_tmp>=power_thresh(i_mice)  & StimOnTime_tmp > 1 & LickEarly_tmp==0 & StimTrials_tmp >0 & Session_type_tmp(:,1)==i_mice & StimDur_tmp>0 & Session_type_tmp(:,4)==SC_side(i_mice));
            elseif i_condition==3
                i_select = (Sample_Delay_tmp==2 & AOM_data_tmp>=power_thresh(i_mice)  & StimOnTime_tmp < 2 & LickEarly_tmp==0 & StimTrials_tmp >0 & Session_type_tmp(:,1)==i_mice & StimDur_tmp>0 & Session_type_tmp(:,4)==SC_side(i_mice));
            elseif i_condition==4
                i_select = (Sample_Delay_tmp==2 & AOM_data_tmp>=power_thresh(i_mice)  & StimOnTime_tmp > 2 & LickEarly_tmp==0 & StimTrials_tmp >0 & Session_type_tmp(:,1)==i_mice & StimDur_tmp>0 & Session_type_tmp(:,4)==SC_side(i_mice));
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

disp('============= perf Contra, p value ==============')
p=sum(delta_perf_btstrp_contra(:,2:5)>=0)/size(delta_perf_btstrp_contra,1)

disp('============= perf ipsi, p value ==============')
p=sum(delta_perf_btstrp_ipsi(:,2:5)<=0)/size(delta_perf_btstrp_ipsi,1)


