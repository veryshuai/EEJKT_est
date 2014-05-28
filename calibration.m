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

pop = [...
%mc best
   0.0968627409766   0.1268472025996   10.7197810687297 0.5080052789763   0.0822185580385   1.2125207260068 5.0343060485555   0.8097175795200   0.0515746775105 0.3370784896322   0.0561120658952   9.5596929422761 0.1977405854859   189.8985217982770;...
   0.1318695142459   0.0947094808226   9.8903670125950 0.6396618909406   0.0793192540062   1.2117586061175 9.6695978460340   0.8083505916956   0.0970437432139 0.3337537574478   0.0739591824742   7.9700195572520 0.3347094431068   182.9484595050208;...
   0.0968627409766   0.1268472025996   10.7197810687297 0.5080052789763   0.0822185580385   1.2125207260068 4.0343060485555   0.8097175795200   0.0515746775105 0.3370784896322   0.0816806661816   9.5596929422761 0.1977405854859   182.4492298503510;...
   %calib best
   0.1366310163926   0.1268472025996   10.5613761750124 0.5080052789763   0.0822185580385   1.2125207260068 1.08592210155339   0.8097175795200   0.0944193209021 0.3370784896322   0.0816806661816   9.0596929422761 0.3439545556451   182.9492298503510;...
   0.1340611886138   0.1302128159986   10.0325890398832 0.5879408294757   0.0622525309197   0.5711740448797 9.2105654192545   0.0420263648961   0.0803515255128 0.2480617991932   0.2461004377215   7.8077669316467 0.3442626654259   166.3152171208176;...
   0.1447534341432   0.1189110004888   8.6129735637322 0.6616151823164   0.0588095440678   0.9936885900071 9.1615265322306   0.2151637864480   0.1763807206260 0.1844713075674   0.2153237303053   7.4273600505579 0.2779843744892   172.7109268948465;...
   0.2521462220765   0.1236548745618   8.6218193102494 0.7289286150099   0.0553133581749   1.5814590436225 9.1912812305071   0.3033580175233   0.0817331611787 0.5836183197545   0.1661335168304   16.2221727715141  0.2750419912764   171.2228381629085;...
   0.2295226979213   0.1383572412773   10.3877590301128 0.6633329881788   0.0614445021089   0.6432478964945 4.9062241524309   0.1927463603387   0.0909765302812 0.3169554240750   0.0578798213703   13.6654637606531 0.2775438664732   177.2906157434897;...
   0.1536413028429   0.1310012440091   9.3277126167340 0.8128799325231   0.0679812303285   0.7788404335537 4.4344885987406   1.1586790703291   0.0819800030894 18.5298845133300  0.3358351041272   170.1710441377340 0.7895076690249;...
%    0.2345468099595   0.1184521509453   10.3584402206035  0.6515722765089  0.0666862159110  0.8653007554242  4.6360100716204  0.4251377177606  0.0602926114751  0.4816971580753  0.2272335228043  17.6733812410265  0.1818407414732  184.5488747841660;...
%    0.2345468099595   0.1184521509453   10.3584402206035  0.6515722765089  0.0666862159110  0.8653007554242  4.6360100716204  0.4251377177606  0.0602926114751  0.4816971580753  0.2272335228043  17.6733812410265  0.1818407414732  184.5488747841660;...
%    0.1083235218887   0.4921435927321   9.4600281850105 0.5103284776595   0.0580380954131   0.5293365525433 8.0844484249265   0.9445379806032   0.3 0.8155180586720 0.3 18.0755054714757   0.108689814316   219.7245451532658;...
 %   0.0580061242556   0.1036625696464   9.6356635240144 0.5729116673072   0.0627897296282   0.9728760370721 8.5736126650738   0.1214101813794   0.4925128567152 0.9440640024788   0.3653144163342   13.6919784473169 0.0861235448786   164.8106946715628;...
 %   0.0580061242556   0.1036625696464   8.6356635240144 0.5729116673072   0.0627897296282   0.9728760370721 8.5736126650738   0.1214101813794   0.4925128567152 0.9440640024788   0.3653144163342   13.6919784473169 0.0861235448786   164.8106946715628;...
 %   0.0580061242556   0.1036625696464   10.6356635240144 0.5729116673072   0.0627897296282   0.9728760370721 8.5736126650738   0.1214101813794   0.4925128567152 0.9440640024788   0.3653144163342   13.6919784473169 0.0861235448786   164.8106946715628;...
   %0.1083235218887   0.4921435927321   9.4600281850105 0.5103284776595   0.0580380954131   0.5293365525433 8.0844484249265   0.9445379806032   0.8155180586720 18.0755054714757   0.108689814316    219.7245451532658;...
   %0.1852718210882   0.2880432861378   9.8829740160325 0.7109566380215   0.0950749999877   0.6777153197839 6.7299645442733   0.3592028580154   1.3022814902186 18.5181492341044   0.3740050586823   253.9086425354636;...
   %0.1852718210882   0.2880432861378   9.8829740160325 0.7109566380215   0.0950749999877   0.6777153197839 6.7299645442733   0.3592028580154   1.3022814902186 18.5181492341044   0.3740050586823   253.9086425354636;...
   %0.1852718210882   0.2880432861378   9.8829740160325 0.7109566380215   0.0950749999877   0.6777153197839 6.7299645442733   0.3592028580154   1.3022814902186 18.5181492341044   0.3740050586823   253.9086425354636;...
   %0.1852718210882   0.2880432861378   9.8829740160325 0.7109566380215   0.0950749999877   0.6777153197839 6.7299645442733   0.3592028580154   1.3022814902186 18.5181492341044   0.3740050586823   253.9086425354636;...
   %0.1852718210882   0.2880432861378   9.8829740160325 0.7109566380215   0.0950749999877   0.6777153197839 6.7299645442733   0.3592028580154   1.3022814902186 18.5181492341044   0.3740050586823   253.9086425354636;...
   %0.1852718210882   0.2880432861378   9.8829740160325 0.7109566380215   0.0950749999877   0.6777153197839 6.7299645442733   0.3592028580154   1.3022814902186 18.5181492341044   0.3740050586823   253.9086425354636;...
   %0.1852107859320   0.2878359811833   9.8829740160325 0.7109566380215   0.0911859161254   0.6777153197839 6.7299645442733   0.3552147708381   1.2996569784998 18.5190654767039  0.3737609180573   253.9080455949462;...
   %0.1852107859320   0.2878359811833   9.8829740160325 0.7109566380215   0.0911859161254   0.6777153197839 6.7299645442733   0.3552147708381   1.2996569784998 18.5190654767039  0.3737609180573   253.9080455949462;...
   %0.1852107859320   0.2878359811833   9.8829740160325 0.7109566380215   0.0911859161254   0.6777153197839 6.7299645442733   0.3552147708381   1.2996569784998 18.5190654767039  0.3737609180573   253.9080455949462;...
   %0.1852107859320   0.2878359811833   9.8829740160325 0.7109566380215   0.0911859161254   0.6777153197839 6.7299645442733   0.3552147708381   1.2996569784998 18.5190654767039  0.3737609180573   253.9080455949462;...
   %0.1852107859320   0.2878359811833   9.8829740160325 0.7109566380215   0.0911859161254   0.6777153197839 6.7299645442733   0.3552147708381   1.2996569784998 18.5190654767039  0.3737609180573   253.9080455949462;...
   %0.1852107859320   0.2878359811833   9.8829740160325 0.7109566380215   0.0911859161254   0.6777153197839 6.7299645442733   0.3552147708381   1.2996569784998 18.5190654767039  0.3737609180573   253.9080455949462;...
        ];

 options=gaoptimset('Display','iter','PopulationSize',12,'Generations',150,...
 'StallTimeLimit',86400,'TimeLimit',3600,'MutationFcn',@mutationadaptfeasible,...
 'FitnessScalingFcn',@fitscalingrank,'InitialPopulation',pop,'UseParallel','always',...
 'PlotFcns',@gaplotbestf,'EliteCount',0);%,'HybridFcn',{@fmincon,fminconoptions});

[X,fval,exitflag,output,population,scores] = ga(@(X) distance(X),14, [],[],[],[],[   0.05;  0.05;    8;    0.5;  .04; 0.5; 1; 0.01; 0.01; 0.01; 0.01; 7; 0.01; 150], [.35;  .4;    11;     0.9;  0.12;    2;    12;    2; 0.5; 2; 0.5; 19; 0.35; 200],[],options); 
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

toc
% matlabpool close
save  calib-10-6-12
diary off
