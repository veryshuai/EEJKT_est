function calc_val(myvalbeta)
% This script calculates the value of the network, using 
% a 30 year window

load('results/val_sim_results'); % get the data

size(sale_f_mat)

dval_init = discount(sum(sale_f_mat(1:15,:),2),myvalbeta);
dval_mature = discount(sum(sale_f_mat(16:30,:),2),myvalbeta);

diff = sum(dval_mature) - sum(dval_init) 

display(['The difference is ', num2str(diff), '.'])

display(sum(sum(sale_f_mat(16:30,:),2)))
display(sum(sum(sale_f_mat(1:15,:),2)))

tot_exp          = log(sum(sale_f_mat,2))

%save('results/val_calc_results')
