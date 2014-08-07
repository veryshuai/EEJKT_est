function val_results = calc_val(myvalbeta)
% This script calculates the value of the network, using 
% a 30 year window

load('results/val_sim_results'); % get the data

val_results = cell(4,1); % this holds the results 

cost_f_mat = cell2mat(cost_f');
net_rev = sum(sale_f_mat,2) - sum(cost_f_mat,2);
tot_rev = sum(sale_f_mat,2)

equil_exp_no = mean(sum(sale_f_mat(26:50,:) > 0,2));

m_net_rev = net_rev / equil_exp_no;
m_tot_rev = tot_rev / equil_exp_no;

dval_init = discount(sum(m_net_rev(1:25,:),2),myvalbeta);
dval_mature = discount(sum(m_net_rev(26:50,:),2),myvalbeta);

diff = sum(dval_mature) - sum(dval_init); 

display(['The difference is ', num2str(diff), ' per active exporter.'])

val_results{1} = diff;
val_results{2} = sale_f_mat;
val_results{3} = cost_f_mat;

%save('results/val_calc_results')
