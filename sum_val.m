% This script takes results of looped network value calculation and makes summary statistics

% Load data
load('results/network_val_calc_results')

% Overall discounted value of the network
diff_list = zeros(size(val_mat));
for k=1:size(val_mat,1)
    diff_list(k) = val_mat{k}{1};
end
m_diff = mean(diff_list);

% Overall value of the network 
diff_list = zeros(size(val_mat));
for k=1:size(val_mat,1)
    diff_list(k) = val_mat{k}{1};
end
m_diff = mean(diff_list);
display(['Discounted value of the network per steady-state active exporter is ' num2str(m_diff) ' 1992 dollars.']); display(' ');

% Costs per year, per searching firm
sale_f = val_mat{1}{2};
cost_f = val_mat{1}{3};
for k=1:size(val_mat,1)
    sale_f = horzcat(sale_f,val_mat{1}{2});
    cost_f = horzcat(cost_f,val_mat{1}{3});
end

% Average number of negative vs positive profit years
ann_prof = sale_f/de - cost_f;
firm_neg = sum(ann_prof < 0, 1);
firm_pos = sum(ann_prof > 0, 1);
neg_year_percent = sum(firm_neg(firm_pos>0)) ./ (sum(firm_neg(firm_pos>0)) + sum(firm_pos(firm_pos>0)));
searchers_with_pos = sum(firm_neg > 0 & firm_pos > 0) / sum(firm_neg > 0);

display([num2str(searchers_with_pos * 100), '% of firms that ever search on the export market ever have positive net profits.']); display(' ');
display(['Firms that have at least one year of positive net profits on the export market have on average ', num2str(sum(firm_neg(firm_pos>0))/sum(firm_pos>0)), ' years with negative net export profit.']); display(' ');
display(['On average, this makes up ', num2str(neg_year_percent), ' of a sometime positive revenue firms lifetime as a potential exporter.']); display(' ');

% Level of search costs
sum_costs = sum(cost_f,1);
display([num2str(sum(sum_costs > 10000) / sum(sum_costs > 0) * 100) '% of firms that ever search spend more than 10,000 USD on search and fixed costs in their lifetime.']); display(' ');

% histograms
h = figure();
annprof_sum = sum(ann_prof);
hist(log(annprof_sum(annprof_sum>0)),100);
title('Lifetime log positive profits')
saveas(h,'results/lifetime_log_pos_profits.jpg')

h = figure();
hist(log(-annprof_sum(annprof_sum<0)),100);
title('Lifetime log negative profits')
saveas(h,'results/lifetime_log_neg_profits.jpg')

h = figure();
hist(log(sum_costs(sum_costs>0)),100);
title('Lifetime log costs')
saveas(h,'results/lifetime_log_costs.jpg')


