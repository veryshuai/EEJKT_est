function dat = read_plot_data(filename)
% reads in data from file, and returns plot data

load(filename);

act_exp          = log(act_exp);
m_spc            = log(m_spc);
m_numclients     = log(m_numclients);
tot_exp          = log(sum(sale_f_mat,2));

dat = {act_exp, m_spc, m_numclients, tot_exp};

