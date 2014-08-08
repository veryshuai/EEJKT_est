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
    sale_f = horzcat(big_list,val_mat{1}{2});
    cost_f = horzcat(big_list,val_mat{1}{3});
end

% Average number of negative vs revenue years
ann_rev = sale_f - cost_f;
firm_neg = sum(ann_rev < 0, 1);
firm_pos = sum(ann_rev > 0, 1);
neg_year_percent = sum(firm_neg(firm_pos>0)) ./ (sum(firm_neg(firm_pos>0)) + sum(firm_pos(firm_pos>0)));
searchers_with_pos = sum(firm_neg > 0 & firm_pos > 0) / sum(firm_neg > 0);

display([num2str(searchers_with_pos * 100), '% of firms that ever search on the export market ever have positive net revenue.']); display(' ');
display(['Firms that have at least one year of positive net revenue on the export market have on average ', num2str(sum(firm_neg(firm_pos>0))/sum(firm_pos>0)), ' years with negative net export revenue.']); display(' ');
display(['On average, this makes up ', num2str(neg_year_percent), ' of a sometime positive revenue firms lifetime as a potential exporter.']); display(' ');

% Level of search costs
sum_costs = sum(cost_f,1);
display([num2str(sum(sum_costs > 100000) / sum(sum_costs > 0) * 100) '% of firms that ever search spend more than 100,000 USD on search and fixed costs in their lifetime.']); display(' ');



