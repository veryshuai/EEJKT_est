%This program calls the fixed parameter settings from SetParams, calls the
%solve function which returns optimal search intensities etc. for every
%state configuration.  Using these intensities, moments are generated in
%the moms function.

%lambda_f is 3 dimensional array with dims [successes,trials,theta_g], each 
%element of the array is a phi_size by x_size matrix giving the optimal
%intensities 

%lambda_h is a 2 dimensional array with dims [theta_g,theta_h], each
%element of the array same as in lambda_f above.

SetParams_noprod;

[lambda_f,lambda_h,pi_tilda_h,pi_tilda_f,c_val_h,c_val_f,punishment] = solve(mm);

[vtran,hazrate,clidist,mstat,mnumex,mavex,mavship,mreg,mexreg,mexshr,mlagereg,mlagdreg,mdeathreg] = moms_nocell(mm,c_val_h,c_val_f,lambda_f,lambda_h);

