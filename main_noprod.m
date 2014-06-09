%This function calls the fixed parameter settings from SetParams, calls the
%solve function which returns optimal search intensities etc. for every
%state configuration.  Using these intensities, moments are generated in
%the moms function.

%lambda_f is 3 dimensional array with dims [successes,trials,theta_g], each 
%element of the array is a phi_size by x_size matrix giving the optimal
%intensities 

%lambda_h is a 2 dimensional array with dims [theta_g,theta_h], each
%element of the array same as in lambda_f above.

% Get parameters
SetParams_noprod;

% Get policy and value functions
[lambda_f,lambda_h,pi_tilda_h,pi_tilda_f,c_val_h,c_val_f,punishment] = solve(mm);

% Simulate model and create loss function statistics
[vtran,hazrate,clidist,mstat,mnumex,mavex,mavship,mreg,mexreg,mexshr,mlagereg,mlagdreg,mdeathreg,simulated_data] = moms_nocell(mm,c_val_h,c_val_f,lambda_f,lambda_h);
