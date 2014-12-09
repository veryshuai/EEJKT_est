%This script uses the results of a single simulation of the model to plot a histogram of profits by measures of productivity, for various exporter ages

%load simulation data
clear; load results/no_learning_sim_results

%format numbers
format long

%Plotting paramters 
c_num = 3; %number of ages to plot (at some point I was calling this cohort, which is confusing)
c_max = 9; %last age to plot

%Create matrix of sales by age 
sale_f_age = ones(18,size(sale_f_mat,2)) * -1; %holds profits by age 
for firm_id=1:size(sale_f_mat,2) %loop through all potential exporters
    nonzero_ind = find(sale_f_mat(:,firm_id)>0,18,'first'); %find indeces of non-zero sales
    if isempty(nonzero_ind) == 0 %only sometime exporters
        export_period = nonzero_ind(end) - nonzero_ind(1) + 1;
        sale_f_age(1:export_period,firm_id) = sale_f_mat(nonzero_ind(1):nonzero_ind(end),firm_id);
    end
end

%Redo the success probability calculation (the one in st_disc is currently screwed up)
succ_prob_correct=zeros(size(sale_f_mat,2),4);
for k=1:size(sale_f_mat,2)
    succ_mat_temp = full(succ_prob_vec{k}); %change the sparse matrix to full, (succ_prob_f,succ_prob_h,succ_ind_f,succ_ind_h)
    if isempty(succ_mat_temp) == 0 %check for emptiness (lots of these...)
        succ_prob_correct(k,:) = succ_mat_temp(end,:); %read from last line into final matrix
    end
end

%Density by cohort
% high productivity 
subplot(2,2,1);
dens_by_cohort(sale_f_age(:,prods>=0.75),c_num,c_max)
title('high productivity');
% low productivity 
subplot(2,2,2);
dens_by_cohort(sale_f_age(:,prods<0.75),c_num,c_max)
title('low productivity')
% high succ
subplot(2,2,3);
mean_succ = mean(succ_prob_correct(succ_prob_correct(:,1)>0,1)); %mean non-zero success probability
dens_by_cohort(sale_f_age(:,succ_prob_correct(:,1)>=mean_succ),c_num,c_max);
title('high success')
% low succ
subplot(2,2,4);
low_suc = sale_f_age(:,succ_prob_correct(:,1)<mean_succ);
dens_by_cohort(low_suc,c_num,c_max);
title('low success')

% Create grid of sales per firm, baseline
sales_grid_cum = zeros(15,7); %results grid, cumulative sales
sales_grid_by_firm = zeros(15,7); %results grid, sales by firm
sale_f_comb = sum(sale_f_mat)'; %combine all years, make it a vector
for succ_ind = 1:7
    for prod_ind = 1:15
        % This sums the sales of all firms
        sales_grid_cum(prod_ind,succ_ind) = sum(sale_f_comb((prods == Phi(prod_ind)) & (succ_prob_correct(:,3) == 8 - succ_ind)));
        % This generates sales by firm
        sales_grid_by_firm(prod_ind,succ_ind) = sales_grid_cum(prod_ind,succ_ind) / max(sum((prods == Phi(prod_ind)) & (succ_prob_correct(:,3) == 8 - succ_ind)),1);
    end
end

sales_grid_cum
sales_grid_by_firm

% Create grid of sales per firm, no learning
sales_grid_cum_nl = zeros(15,7); %results grid, cumulative sales
sales_grid_by_firm_nl = zeros(15,7); %results grid, sales by firm
sale_f_comb = sum(sale_h_mat)'; %combine all years, make it a vector
for succ_ind = 1:7
    for prod_ind = 1:15
        % This sums the sales of all firms
        sales_grid_cum_nl(prod_ind,succ_ind) = sum(sale_f_comb((prods == Phi(prod_ind)) & (succ_prob_correct(:,4) == 8 - succ_ind)));
        % This generates sales by firm
        sales_grid_by_firm_nl(prod_ind,succ_ind) = sales_grid_cum_nl(prod_ind,succ_ind) / max(sum((prods == Phi(prod_ind)) & (succ_prob_correct(:,4) == 8 - succ_ind)),1);
    end
end

sales_grid_cum_nl
sales_grid_by_firm_nl
