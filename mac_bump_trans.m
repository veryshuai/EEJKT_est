function out = mac_bump_trans(garbage)
out = garbage;
% This script runs the counterfactual of increasing the foreign macro process

increase = 1.2; %factor with which to multiply foreign macro shocks

%estimated parameters
X = [...
   0.0338022054957   0.2671671430852   11.3440817512744 0.5121388785243   0.0870454860595   0.7155086649572 3.1608659441438   0.5320902587029   0.0868923804355 8.8391588460077   0.2983646166548   110.3891766498293 0.6479607484458;...
        ];

rng(80085);

tic

format long;

%read in parameter names from the X vector
scale_h    =  X(3);
ag         =  .5;
bg         =  .5;
lnF        =  scale_h + log(X(1));
delta      =  X(2);
scale_f_orig    =  scale_h + log(X(4));
scale_f_new    =  scale_h + log(X(4) * increase);
beta       =  X(5);
ah         =  X(6);
bh         =  X(7);
L_p        =  0; %zero would require changing the value function solver...
D_p        =  0;
L_z        =  X(8);
D_z        =  X(9);
L_b        =  X(10);
alp        =  0;
gam        =  X(11)*(1+beta)/beta;
cs         =  X(12);
sig_p      =  X(13);

main_trans;  %this is the main program (to run simulation, generate moments, etc)

%% Targets

vtranData     = [0.618; 0.321; 0.048; 0.013; 0.271; 0.375; 0.241; 0.113]; %Taken from CJ_cohorts_6-15-2012.pdf in S://moments
mhazData      = [0.694; 0.515; 0.450; 0.424;.389]; %Taken from CJ_cohorts_6-15-2012.pdf in S://moments
clidistData   = [-1.667; -0.097;0.0663];  %Taken from file "Copy of client distribution regression.xslx" in S://moments, ultimately from CJ's data, not sure about the date of email
%mstatData     = [9.746;2.107]; %Taken from CJ_cohorts_6-15-2012.pdf in S://moments, and calculated from CJ's email of 12/18/2012 respectively
mnumexData    = [0.709; 0.383; 0.300; 0.263; 0.293]; %Taken from CJ_cohorts_6-15-2012.pdf in S://moments
mavexData     = [8.960; 10.018; 10.231; 10.370; 10.473]; % From CJ's email of March 29, 2013 
mavshipData   = 4.824; % dropping top 5%, average number of shipments per year per client via CJ's email of 6-9-2012
mregData      = [2.677;-0.143;2.180]; %ncoefficients from reg of log(sales per client) on log no. of clients, (log no. of clients)^2, and root MSE, from file "EINS as IDs results2.docx" in CJ's email of 6-11-2012
mexregData    = [0.727;2.167]; % log exports on log domestic sales, coefficient and lMSE from eam_moms_out( DIAN only 8 Feb 2013).log file in Marcela's email of 2-11-2013, on S: drive in 'Marcela Stuff'.
mexshrData    = 0.299; % share of exporting firms, from eam_moms_out( DIAN only 8 Feb 2013).log file in Marcela's email of 2-11-2013, on S: drive in 'Marcela Stuff'.
mlageregData  = [0.811;0.233;1.124]; %This is now the match level ar1 from CJ's email of 3-27-2013, log match sales on lagged match sales, first year, and MSE, last model with HS dummies
mlagdregData  = [0.976;0.462]; %regression coefficient and MSE of log domestic sales on log lagged domestic sales in Marcela's email of 3-22-2013 
mdeathregData = [1.174;0.166;-0.070;0.453]; %intercept, coefficient on first yr, and ln sales, dependent death dummy, source is last model (annual with HS dummies) in CJ's email of 3-27-2013.

%% Weights

%standard errors (see moments folder for details of calculations)
vtranSE     = [0.0025; 0.0024; 0.0011; 0.0006; 0.0071; 0.0066; 0.0059; 0.0043]; %Calculated from CJ_cohorts_6-15-2012.pdf in S://moments
mhazSE      = [0.00230; 0.00484; 0.00674; 0.00865; 0.01066]; %Calculated from CJ_cohorts_6-15-2012.pdf in S://moments
clidistSE   = [0.0566;0.0132;sqrt(2)*0.0663^2/sqrt(36)]; %Calculated from file "Copy of client distribution regression.xslx" in S://moments
%mstatSE     = [2.107/sqrt(48133);sqrt(2)*2.197^2/sqrt(48133)]; %all calculated from CJ's email of 12-18-2012.
mnumexSE    = [0.00287; 0.00569; 0.00683; 0.00784;0.00861]; %Calculated from CJ_cohorts_6-15-2012.pdf in S://moments
mavexSE     = [1.821/sqrt(22838); 1.974/sqrt(6484); 2.044/sqrt(3922); 2.063/sqrt(2713);2.063/sqrt(1953)]; %from CJ's email of 3-27-2013.  The two 2.063's are intentional, with rounding they are the same in the data.
mavshipSE   = 7.142/sqrt(67432); % via CJ's email of 6-7-2012, total number of pairs is from CJ's email of 5-17-2012.
mregSE      = [0.15195;0.09592;sqrt(2)*2.180^2/sqrt(6618)];  % Calculated from from file "EINS as IDs results2.docx" in CJ's email of 6-11-2012
mexregSE    = [0.0065;sqrt(2)*2.127^2/sqrt(34330)]; % calculated from eam_moms_out.log file in Marcela's email of 9-24-2012
mexshrSE    = 0.466/sqrt(115334); % calculated from eam_moms_out.log file in Marcela's email of 9-24-2012
%mlagsregSE  = 0.00367/sqrt(20223); % calculated from in-house transactions data, .do file/log is in S: moments folder.
mlageregSE  = [0.00347;0.0143;sqrt(2)*1.124^2/sqrt(33133)]; %from CJ's 3-27-2013 email, last model
mlagdregSE  = [0.000862;sqrt(2)*0.462^2/sqrt(99302)]; %regression coefficient and MSE of log dom sales on log lagged dom sales, in-house 80's colombia 
mdeathregSE = [0.0265;0.00348;0.000847;sqrt(2)*0.453^2/sqrt(81383)]; %CJ's email of 3-27-2013

Data = cat(1,vtranData,mhazData,clidistData,mnumexData,mavshipData,mavexData,mregData,mexregData,mexshrData,mlageregData,mlagdregData,mdeathregData);

W = zeros(size(Data,1));
W(1:size(Data,1)+1:end) = (1./cat(1,vtranSE,mhazSE,clidistSE,mnumexSE,mavshipSE,mavexSE,mregSE,mexregSE,mexshrSE,mlageregSE,mlagdregSE,mdeathregSE)).^2;

%% Realizations

% Model = cat(1,vtran,hazrate,clidist,mstat,mnumex,mavship,mavex,mreg,mexreg,mexshr);
Model = cat(1,vtran,hazrate,clidist,mnumex,mavship,mavex,mreg,mexreg,mexshr,mlagereg,mlagdreg,mdeathreg);

try
    error = Data-Model;

    D = error'*W*error;

    nanflag = isnan(D); 
    if nanflag>0;
        D = D * 10; 
    end

    punishment
    D = D*(1+punishment);
    
    mmm = cat(2,Data,Model)
    X %print out current parameter guess
    D

    %Old D report
    Old_D = norm(error)/norm(Data)

    last_loop_run_time = toc/60
catch err 
    D = 1e12;
end

dlmwrite('mac_bump.txt', Model)
dlmwrite('mac_bump.txt', D, '-append')
dlmwrite('mac_bump.txt', Old_D, '-append')

end
