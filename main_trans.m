%This program calls the fixed parameter settings from SetParams, calls the
%solve function which returns optimal search intensities etc. for every
%state configuration.  Using these intensities, moments are generated in
%the moms function.

%lambda_f is 3 dimensional array with dims [successes,trials,theta_g], each 
%element of the array is a phi_size by x_size matrix giving the optimal
%intensities 

%lambda_h is a 2 dimensional array with dims [theta_g,theta_h], each
%element of the array same as in lambda_f above.

scale_f = scale_f_orig
SetParams_trans;
[lambda_f_orig,lambda_h_orig,pi_tilda_h,pi_tilda_f,c_val_h_orig,c_val_f_orig,punishment] = solve(mm);

scale_f = scale_f_new
SetParams_trans;
[lambda_f_new,lambda_h_new,pi_tilda_h,pi_tilda_f,c_val_h_new,c_val_f_new,punishment] = solve(mm);

save policy_trans

[vtran,hazrate,clidist,mstat,mnumex,mavex,mavship,mreg,mexreg,mexshr,mlagereg,mlagdreg,mdeathreg] = moms_trans(mm,c_val_h_orig,c_val_f_orig,c_val_h_new,c_val_f_new,lambda_f_orig,lambda_h_orig,lambda_f_new,lambda_h_new,scale_f_orig,scale_f_new);

