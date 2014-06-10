% This script returns simulated data
% for use in creating plots and counterfactuals

% Parallel setup
clc
if matlabpool('size')~=12 %if pool not equal to 12, open 12
   if matlabpool('size')>0
     matlabpool close
   end
   matlabpool open 12
end

% Add correct root folder
addpath(genpath('/gpfs/home/dcj138/work/Colombia-Project/Colombia-Project/'));

% Put numbers in long format for printing
format long;

% random seed
rng(80085);

% Initial condition
X = [...
   0.0338022054957   0.2671671430852   11.3440817512744 0.5121388785243   0.0870454860595   0.7155086649572 3.1608659441438   0.5320902587029   0.0868923804355 8.8391588460077   0.2983646166548   110.3891766498293 0.6479607484458;...
        ];

[D,W,error,simulated_data] = distance_noprod(X)

% lnF         =  scale_h+log(X(1));
% delta       =  X(2);
% scale_h     =  X(3);
% scale_f     =  scale_h + log(X(4));
% beta        =  X(5);
% a           =  X(6);
% b           =  X(7);
% L_p         =  X(8);
% D_p         =  X(9);
% L_z         =  X(10);
% D_z         =  X(11);
% L_b         =  X(12);
% gam         =  X(13)/beta;
% cost scalar =  X(14);

% End timer
toc

% Close parallel
matlabpool close

% Save results
save est_results;
