function D = distance(X,shk)
% This function takes a set of parameters and returns the distance between
% the moments in the data and the moments generated by the parameter set.

%last updated 2011/8/23 by David Jinkins

tic

%read in parameter names from the X vector
scale_h    =  X(5);
ag        =  X(1);
bg        =  X(2);
lnF        =  scale_h + log(X(3));
delta      =  X(4);
scale_f    =  scale_h + log(X(6));
beta        =  X(7);
ah         =  X(8);
bh         =  X(9);
L_p        =  X(10);
L_z        =  X(11);
L_b        =  X(12);
alp        =  X(13);

%display(X);

main1;  %this is the main program (to run simulation, generate moments, etc)

%% Targets


vtranData     = [0.614; 0.332; 0.043; 0.011; 0.262; 0.375; 0.249; 0.115]; %This is chance of moving [1->0;1->1;1->2;1->>3,2->0;2->1;2->2;2->>3]
mhazData      = [0.713; 0.378; 0.289; 0.281; 0.217]; %Probability of match dying conditional on surviving 1,2,3,4 and 5 years
clidistData   = [.843; .100; .057];  %percentage of firms with 1,2,3,4, and 5 clients
mstatData     = [0.1870; 0.366; 0.4947; 0.3045]; %aggregate stats, see moms.m for specifics
mnumexData    = [0.507; 0.441; 0.422; 0.390]; % 2-5 yr old cohort membership relative to one year olds year
mavexData     = [0.312; 0.237; 0.217; 0.191]; % 2-5 yr. old inverse sales per exporter x first year sales per exporter
mavshipData   = [1/8;1/8;4/8;14/100;53/100]; %[10%,25%,50%,75%,90%] percentile of number of shipments per year divided by 8 and 100 
mregData      = [.558;-0.00589*(-100);3.493*(1/5)]; %normalized coefficients from reg of log(sales per client) on no. of clients, no. of clients^2, and MSE, from file reg2 lsales = #buyers.docx in cj_table_and_moment_data folder

Data = cat(1,vtranData,mhazData,clidistData,mstatData,mnumexData,mavexData,mavshipData,mregData);
%% Realizations

Model = cat(1,vtran,hazrate,clidist,mstat,mnumex,mavex,mavship,mreg);
   
    D = norm(Data-Model)/norm(Data);
    D_report = D;
    
    nanflag = isnan(D); 
    if nanflag>0;
        D = 100; 
    end
    
    display('vtran(8),mhaz(5),clidist(3),mstat(4),mnumex(4),mavex(4),mavship(5),mreg(3)'); 
    mmm = cat(2,Data,Model)
    X %print out current parameter guess
    D_report


toc
end
