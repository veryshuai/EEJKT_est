%This script uses the results of a single simulation of the model to plot a histogram of profits by measures of productivity, for various cohort years

%load simulation data
clear; load results/cf_sim_results

%Create matrix of sales by cohort
sale_f_coh = ones(18,size(sale_f_mat,2)) * -1; %holds profits by cohorts
for firm_id=1:size(sale_f_mat,2) %loop through all potential exporters
    nonzero_ind = find(sale_f_mat(:,firm_id)>0,18,'first'); %find indeces of non-zero sales
    if isempty(nonzero_ind) == 0 %only sometime exporters
        export_period = nonzero_ind(end) - nonzero_ind(1) + 1;
        sale_f_coh(1:export_period,firm_id) = sale_f_mat(nonzero_ind(1):nonzero_ind(end),firm_id);
    end
end

