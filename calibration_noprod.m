% Calibration
%This script initiates the genetic algorithm to solve for
%the parameters of the model.  Initial parameter guesses can be set here,
%as well as optional settings for the genetic algorithm.

% Parallel setup
clc
if matlabpool('size')~=12 %if pool not equal to 12, open 12
   if matlabpool('size')>0
     matlabpool close
   end
   matlabpool open 12
end

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
   0.1943930477380   0.2453104276277   9.0369335893045 0.2662888213532   0.0416788862843   0.6809826297676 4.4510131294695   0.7066529039047   0.0934049031329 1.2680600397205   0.3686143534469   398.4058091796875 0.3925683743937;...
   0.1943930477380   0.2453104276277   9.0369335893045 0.2662888213532   0.0416788862843   0.6809826297676 4.4510131294695   0.7066529039047   0.0934049031329 1.2680600397205   0.3686143534469   398.4058091796875 0.3925683743937;...
   0.1943930477380   0.2453104276277   9.0369335893045 0.2662888213532   0.0416788862843   0.6809826297676 4.4510131294695   0.7066529039047   0.0934049031329 1.2680600397205   0.3686143534469   398.4058091796875 0.3925683743937;...
   0.1943930477380   0.2453104276277   9.0369335893045 0.2662888213532   0.0416788862843   0.6809826297676 4.4510131294695   0.7066529039047   0.0934049031329 1.2680600397205   0.3686143534469   398.4058091796875 0.3925683743937;...
   0.1943930477380   0.2453104276277   9.0369335893045 0.2662888213532   0.0416788862843   0.6809826297676 4.4510131294695   0.7066529039047   0.0934049031329 1.2680600397205   0.3686143534469   398.4058091796875 0.3925683743937;...
   0.1943930477380   0.2453104276277   9.0369335893045 0.2662888213532   0.0416788862843   0.6809826297676 4.4510131294695   0.7066529039047   0.0934049031329 1.2680600397205   0.3686143534469   398.4058091796875 0.3925683743937;...
    ];

% Options for genetic algorithm
options = gaoptimset('Display','iter','PopulationSize',12,'Generations',150,...
'StallTimeLimit',86400,'TimeLimit',3600,'MutationFcn',@mutationadaptfeasible,...
'FitnessScalingFcn',@fitscalingrank,'InitialPopulation',pop,'UseParallel','always',...
'PlotFcns',@gaplotbestf,'EliteCount',0);%,'HybridFcn',{@fmincon,fminconoptions});

% Call estimation routine
[X,fval,exitflag,output,population,scores] = ga(@(X) distance_noprod(X),13, [],[],[],[],[0.1; 0.1; 8; 0.1; 0.01; 0.01; 0.01;  0.03; 0.01; 1; 0; 200; .1], [.5;  0.5; 11.0; 0.5; 0.3; 5; 5; 1; 0.4; 5; 0.5; 550; 1],[],options);  
%[X,fval,exitflag,output,population,scores] = ga(@(X) distance_noprod(X),13, [],[],[],[],[   0.005;  0.01;    7;    0.1;  .01; 0.1; 0.1;  0.01; 0.01; 1; 0.00; 100; .01], [.5;  1;    12;     1;  0.2;    3;  13; 2; 1; 15; 0.4; 550; 2],[],options);  
   
% lnF         =  scale_h+log(X(1));
% delta       =  X(2);
% scale_h     =  X(3);
% scale_f     =  scale_h + log(X(4));
% beta        =  X(5);
% a           =  X(6);
% b           =  X(7);
% L_z         =  X(8);
% D_z         =  X(9);
% L_b         =  X(10);
% gam         =  X(11)*(1+beta)/beta;
% cost scalar =  X(12);
% sig p       =  X(13);

% End timer
toc

% Close parallel
matlabpool close

% Save results
save est_results;
