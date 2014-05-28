%SetParams


th_ind              = zeros(2000,3);
mu_f                = zeros(2000,1);
mu_h                = zeros(2000,1);
sp_p                = zeros(2000,200);
lambda_f            = zeros(21,21,3,31,15,15);
lambda_h            = zeros(3,7,31,15,15);
c_val_h             = zeros(15,15,15);
c_val_f             = zeros(15,15,15);
burn                = zeros(1,1);
delta               = zeros(1,1);
d                   = zeros(1,1);
S                   = zeros(1,1);
n_size              = zeros(1,1);
net_size            = zeros(1,1);
Z                   = zeros(15,1);
Phi                 = zeros(15,1);
X_f                 = zeros(15,1);
X_h                 = zeros(15,1);
actual_h            = zeros(18,2);
actual_f            = zeros(20,2);
L_b                 = zeros(1,1);
L_z                 = zeros(1,1);
L_f                 = zeros(1,1);
L_h                 = zeros(1,1);
erg_pz              = zeros(15,1);
erg_pp              = zeros(15,1);
maxc                = zeros(1,1);
max_client_prod     = zeros(1,1);
mult_match_max      = zeros(1,1);
mms                 = zeros(1,1);
j                   = zeros(1,1);
max_mat_violation   = zeros(1,1);
match_violation     = zeros(1,1);
violation           = zeros(1,1);
match_number_violation= zeros(1,1);
no_more_rands       = zeros(1,1);
breakflag           = zeros(1,1);
x_size              = zeros(1,1);
Phi_size            = zeros(1,1);
z_size              = zeros(1,1);
TT                  = zeros(1,1);
cum_erg_pz          = zeros(15,1);
cum_erg_pp          = zeros(15,1);
cum_sp_p            = zeros(2000,200);

singlefirm(j,max_mat_violation,match_violation,violation,match_number_violation,no_more_rands,breakflag,x_size,Phi_size,z_size,TT,cum_erg_pz,cum_erg_pp,cum_sp_p,th_ind,mu_h,mu_f,sp_p,lambda_f,lambda_h,c_val_h,c_val_f,burn,delta,d,S,n_size,net_size,Z,Phi,X_f,X_h,actual_h,actual_f,L_b,L_z,L_f,L_h,erg_pz,erg_pp,maxc,max_client_prod,mult_match_max,mms);
