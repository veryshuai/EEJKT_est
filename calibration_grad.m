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

% format long
format long;

% random seed, needed to call pop values
rng(80085);
myStream = rand(10000000,1);

pop = [...
   0.1847434855170   0.2878359811833   9.8822349184373 0.7109566380215   0.0909589416381   0.7173997437809 6.7299645442733   0.3545653186286   1.2989560278772 19.1326667984615   0.3734662326935   254.9079988649047;...
        ];

 options=optimset('Display','iter');

[X,fval] = fmincon(@(X) distance(X,myStream),pop, [],[],[],[],[0.05;  0.1;    8;    0.5;  .05; 0.5; 1; 0.01; 0.01;  13; 0.3; 180], [.3;  .5;    12;     0.9;  0.10;    2;    12;    2;  2; 25; 0.70; 260],[],options);

% lnF         =  scale_h+log(X(1));
% delta       =  X(2);
% scale_h     =  X(3);
% scale_f     =  scale_h + log(X(4));
% beta        =  X(5);
% bmh         =  X(6);
% bvh         =  X(7) * X(6) * (1-X(6));
% L_p         =  X(8);
% L_z         =  X(9);
% L_b         =  X(10);
% gam         =  X(11)/beta;
% cost scalar =  X(12);

toc
% matlabpool close
save  calib-10-6-12
diary off
