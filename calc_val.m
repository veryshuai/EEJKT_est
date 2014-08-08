function val_results = calc_val(myvalbeta,demand_elas)
% This script calculates the value of the network, using 
% a 30 year window

load('results/val_sim_results'); % get the data

val_results = cell(4,1); % this holds the results 

cost_f_mat = cell2mat(cost_f');
net_prof = sum(sale_f_mat,2)/demand_elas - sum(cost_f_mat,2);
tot_prof = sum(sale_f_mat,2)/demand_elas

equil_exp_no = mean(sum(sale_f_mat(30:59,:) > 0,2));

m_net_prof = net_prof / equil_exp_no;
m_tot_prof = tot_prof / equil_exp_no;

dval_init = discount(sum(m_net_prof(1:29,:),2),myvalbeta);
dval_mature = discount(sum(m_net_prof(30:59,:),2),myvalbeta);

diff = sum(dval_mature) - sum(dval_init); 

display(['The difference is ', num2str(diff), ' per active exporter.'])

val_results{1} = diff;
val_results{2} = sale_f_mat;
val_results{3} = cost_f_mat;

%save('results/val_calc_results')
