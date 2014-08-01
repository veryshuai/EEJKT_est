function calc_val(b)
% This script calculates the value of the network, using 
% a 30 year window

load('results/val_sim_results') % get the data

dval_init = discount(sum(sale_f_mat(1:30,:),2),b)
dval_mature = discount(sum(sale_f_mat(31:60,:),2),b)

diff = sum(dval_mature) - sum(dval_init) 

save('results/val_calc_results')
