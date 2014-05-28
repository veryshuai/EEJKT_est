% This script converts the estimates from the calibration routine into the parameters in the model, including SEs

load se-5-3-2013

%lnF
lnF     = X(3) + log(X(1));
lnF_d1  = 1/X(1);
lnF_d3  = 1;
lnF_mid = [varcov(1,1),varcov(1,3);varcov(1,3),varcov(3,3)];
lnF_var = [lnF_d1 lnF_d3] * lnF_mid * [lnF_d1 lnF_d3]';
lnF_se  = sqrt(lnF_var) 
display('Log fixed costs');
display(lnF);
display(lnF_se);
display(' ');

%scale_f
scale_f     = X(3) + log(X(4));
scale_f_d3  = 1;
scale_f_d4  = 1/X(4);
scale_f_mid = [varcov(3,3),varcov(3,4);varcov(3,4),varcov(4,4)];
scale_f_var = [scale_f_d3 scale_f_d4] * scale_f_mid * [scale_f_d3 scale_f_d4]';
scale_f_se  = sqrt(scale_f_var) 
display('foreign scale');
display(scale_f);
display(scale_f_se);
display(' ');

%gamma
gamma     = X(11) * (1+X(5))/X(5);
gamma_d11 = (1+X(5))/X(5);
gamma_d5  = -X(11)/(X(5))^2;
gamma_mid = [varcov(11,11),varcov(5,11);varcov(5,11),varcov(5,5)];
gamma_var = [gamma_d11 gamma_d5] * gamma_mid * [gamma_d11 gamma_d5]';
gamma_se  = sqrt(scale_f_var) 
display('gamma');
display(gamma);
display(gamma_se);
display(' ');

