% This script takes results of looped network value calculation and makes summary statistics

% Load data
load('results/network_val_calc_results')

% Overall discounted value of the network
diff_list = zeros(size(val_mat));
for k=1:size(val_mat,1)
    diff_list(k) = val_mat{k}{1};
end
m_diff = mean(diff_list);





