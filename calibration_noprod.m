% Calibration
%This script initiates the genetic algorithm to solve for
%the parameters of the model.  Initial parameter guesses can be set here,
%as well as optional settings for the genetic algorithm.

% Parallel setup
% clc
% if matlabpool('size')~=12 %if pool not equal to 12, open 12
%    if matlabpool('size')>0
%      matlabpool close
%    end
%    matlabpool open 12
% end

% Start timer
tic;

% Add correct root folder
addpath(genpath('/gpfs/home/dcj138/work/Colombia-Project/Colombia-Project/'));

% Put numbers in long format for printing
format long;

% random seed
rng(80085);

% Initial condition
pop = [...
   0.0338022054957   0.2671671430852   11.3440817512744 0.5121388785243   0.0870454860595   0.7155086649572 3.1608659441438   0.5320902587029   0.0868923804355 8.8357713948358   0.2983646166548   111.4990475604738 0.6496979258178;...
   0.0338022054957   0.2671671430852   11.3440817512744 0.5121388785243   0.0870454860595   0.7155086649572 3.1608659441438   0.5320902587029   0.0868923804355 8.8357713948358   0.2983646166548   111.4990475604738 0.6496979258178;...
   0.0338022054957   0.2671671430852   11.3440817512744 0.5121388785243   0.0870454860595   0.7155086649572 3.1608659441438   0.5320902587029   0.0868923804355 8.8357713948358   0.2983646166548   111.4990475604738 0.6496979258178;...
   0.0338022054957   0.2671671430852   11.3440817512744 0.5121388785243   0.0870454860595   0.7155086649572 3.1608659441438   0.5320902587029   0.0868923804355 8.8357713948358   0.2983646166548   111.4990475604738 0.6496979258178;...
   0.0797622055857   0.1199762182097   11.6147635503327 0.3180748359041   0.0568771088557   0.7870806852870 0.9626329038517   0.4075632811638   0.0868923804355 9.8391588460077                   0   174.9469746301988 0.8351230612318;...
   0.0797622055857   0.1199762182097   12.6147635503327 0.3180748359041   0.0568771088557   0.7870806852870 0.9626329038517   0.7825632811638   0.0868923804355 10.0266588460077                   0   174.9469746301988 0.9601230612318;...
        ];

% Options for genetic algorithm
options = gaoptimset('Display','iter','PopulationSize',12,'Generations',150,...
'StallTimeLimit',86400,'TimeLimit',3600,'MutationFcn',@mutationadaptfeasible,...
'FitnessScalingFcn',@fitscalingrank,'InitialPopulation',pop,'UseParallel','always',...
'PlotFcns',@gaplotbestf,'EliteCount',0);%,'HybridFcn',{@fmincon,fminconoptions});

% Call estimation routine
[X,fval,exitflag,output,population,scores] = ga(@(X) distance_noprod(X),13, [],[],[],[],[   0.005;  0.01;    9;    0.1;  .02; 0.1; 0.1;  0.01; 0.01; 1; 0.00; 30; .01], [.5;  1;    15;     1;  0.10;    3;  15; 2; 1; 15; 0; 250; 2],[],options);  % gamma bounds were [.01 .6]
   
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
save est_results, replace
