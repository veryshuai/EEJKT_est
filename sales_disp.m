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

%Redo the success probability calculation (the one in st_disc is currently screwed up)
succ_prob_correct=zeros(size(sale_f_mat,2),4);
for k=1:size(sale_f_mat,2)
    succ_mat_temp = full(succ_prob_vec{k}); %change the sparse matrix to full
    if isempty(succ_mat_temp) == 0 %check for emptiness (lots of these...)
        succ_prob_correct(k,:) = succ_mat_temp(end,:); %read from last line into final matrix
    end
end

%Density by cohort
% high productivity 
subplot(2,2,1);
dens_by_cohort(sale_f_coh(:,prods>=2))
% low productivity 
subplot(2,2,2);
dens_by_cohort(sale_f_coh(:,prods<2))
% high succ
subplot(2,2,3);
mean_succ = mean(succ_prob_correct(succ_prob_correct(:,1)>0,1)); %mean non-zero success probability
dens_by_cohort(sale_f_coh(:,succ_prob_correct(:,1)>=mean_succ));
% low succ
subplot(2,2,4);
dens_by_cohort(sale_f_coh(:,succ_prob_correct(:,1)<mean_succ));

