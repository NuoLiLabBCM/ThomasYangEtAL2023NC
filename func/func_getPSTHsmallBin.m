function [PSTH time] = func_getPSTHsmallBin(SpikeTimes, PSTH_StartTime, PSTH_EndTime)

% 
% SpikeTimes -- {n_rep,1}
% 

if nargin == 1
    PSTH_StartTime = -.52;
    PSTH_EndTime = 5.020;
end

time = PSTH_StartTime-.05:.05:PSTH_EndTime+.05;


n_rep = size(SpikeTimes,1);
total_counts = 0;
for i_rep = 1:n_rep
    
    [counts] = hist(SpikeTimes{i_rep,1},time);
    total_counts = total_counts+counts/n_rep;
    
end


PSTH = total_counts/.05;


time = time(2:end-2);
PSTH = PSTH(2:end-2);

return