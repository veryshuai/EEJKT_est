% Given a Monte Carlo simulation of firms done in the bootstrap.m script (using the saved file bcov.mat), this script re-estimates the model using generated moments.

clc

load bcov.mat
fakemoments = moments(:,1);
clearvars -except fakemoments;

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
   0.1301176595241   0.5399979003269   11.6022592260359 0.9816148700071   0.0794507619431   0.8964487983315 1.6207439301099   0.9212152597775   0.0471184375495 13.9265096878654   0.3371862600220   226.7746985032681 0.7404109597207;...
   0.1301176595241   0.5399979003269   11.6022592260359 0.9816148700071   0.0794507619431   0.8964487983315 1.6207439301099   0.9212152597775   0.0471184375495 13.9265096878654   0.3371862600220   226.7746985032681 0.7404109597207;...
   0.1341299908906   0.4102899572148   11.1376458090273 0.9867625740705   0.0811061439927   1.8490704345995 2.6168912631247   1.7775849727227   0.0790358402786 18.3721949798804   0.3103330808340   226.2880726083640 1.0439698055737;...
      ];

options=gaoptimset('Display','iter','PopulationSize',12,'Generations',150,...
'StallTimeLimit',86400,'TimeLimit',3600,'MutationFcn',@mutationadaptfeasible,...
'FitnessScalingFcn',@fitscalingrank,'InitialPopulation',pop,'UseParallel','always',...
 'PlotFcns',@gaplotbestf,'EliteCount',0);%,'HybridFcn',{@fmincon,fminconoptions});

[X,fval,exitflag,output,population,scores] = ga(@(X) distance_mc(X,fakemoments),13, [],[],[],[],[   0.05;  0.05;    9;    0.5;  .04; 0.5; 1;  0.01; 0.01; 13; 0.01; 140; .01], [.35;  .6;    12;     1;  0.10;    3;  8; 2; 0.2; 24; .6; 250; 2],[],options); 

% bmg         =  X(1);
% bsdh        =  X(2);
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
% alp         =  X(13);
% gam         =  X(14)/beta;
% cost scalar =  X(15);

toc
matlabpool close
save  montecarlo-12-9-12


