% This script calls the policy plot function

clear 
load policy_trans_orig
policy_plot('original',lambda_f,0,lambda_f)

clear
load policy_trans_no_net
policy_plot('no network',lambda_f_new,0,lambda_f_new)

clear
load policy_trans
policy_plot('macro shock increase',lambda_f_new,0,lambda_f_new)
policy_plot('macro shock increase',lambda_f_orig,1,lambda_f_new)

clear
load policy_trans_cost_dec
policy_plot('fixed cost decrease',lambda_f_new,0,lambda_f_new)
policy_plot('fixed cost decrease',lambda_f_orig,1,lambda_f_new)

clear
load policy_trans_search_dec
policy_plot('search cost decrease',lambda_f_new,0,lambda_f_new)
policy_plot('search cost decrease',lambda_f_orig,1,lambda_f_new)

clear
load policy_trans_red_var
policy_plot('macro shock variance decrease',lambda_f_new,0,lambda_f_new)
policy_plot('macro shock variance decrease',lambda_f_orig,1,lambda_f_new)
