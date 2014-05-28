% Calibration
%This script initiates the genetic algorithm to solve for
%the parameters of the model.  Initial parameter guesses can be set here,
%as well as optional settings for the genetic algorithm.

clc
if matlabpool('size')~=12
  if matlabpool('size')>0
    matlabpool close
  end
  matlabpool open 12
end

tic;
addpath(genpath('/gpfs/home/dcj138/work/Colombia-Project/Colombia-Project/'));

rng(3707);

pop = [...
    % 0.1525    0.1512   10.8005    0.6419    0.0318    0.1862    0.8258 0.8195    0.1455   12.5996    0.1569   92.1079;...
    0.1067    0.1352   10.3982    0.4586    0.0469    0.0970 0.7038    0.9592    0.1405   13.4647    0.0768   91.3634;...
];

% options=gaoptimset('Display','iter','PopulationSize',13,'Generations',150,...
% 'StallTimeLimit',86400,'TimeLimit',3600,'MutationFcn',@mutationadaptfeasible,...
% 'FitnessScalingFcn',@fitscalingrank,'InitialPopulation',pop,'UseParallel','always',...
% 'PlotFcns',@gaplotbestf,'EliteCount',0);%,'HybridFcn',{@fmincon,fminconoptions});

options=saoptimset('Display','iter','PlotFcns',@saplotbestf);

%[X,fval,exitflag,output,population,scores] = ga(@(X) distance(X),12, [],[],[],[],[   0;  .1;    10;    .3;  .03; 0.01; 0.01; 0.01;  0.01;  5; 0; 50], [.5;  .5;    16;     1;  .07;    .95;    3;    3;     3; 16; .4; 200],[],options);
[X,fval,exitflag,output] = simulannealbnd(@(X) distance_sa(X),pop, [   0;  .1;    10;    .3;  .04; 0.01; 0.01; 0.01;  0.01;  9; 0; 50], [.5;  .7;    16;     .9;  .09;    .99;    3;    3;     3; 16; .4; 220],options);

% lnF         =  scale_h+log(X(3));
% delta       =  X(4);
% scale_h     =  X(5);
% scale_f     =  scale_h + log(X(6));
% beta        =  X(7);
% bmh         =  X(8);
% bsdh        =  X(9);
% L_p         =  X(10);
% L_z         =  X(11);
% L_b         =  X(12);
% gam         =  X(14)/beta;
% cost scalar =  X(15);

toc
% matlabpool close
save  sa-9-9-12
diary off
