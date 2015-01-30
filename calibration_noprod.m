function [] = calibration_noprod(pop, varargin)
%This function initiates the genetic algorithm to solve for
%the parameters of the model.  Initial parameter guesses can be set here,
%as well as optional settings for the genetic algorithm.

    % set defaults for optional inputs
    optargs = {0};

    % overwrite defaults with user input
    numvarargs = size(varargin, 2);
    optargs(1:numvarargs) = varargin;

    % memorable variable names
    rst = optargs{1}; % parameter restriction?

    % what type of parameter restriction?
    rst_type = 'non';
    if rst == 1
        rst_type = query(['Which type of parameter restriction?' char(10) 'nln (no learning)' char(10) 'nnt (no network effect)?'],{'nln','nnt'});
    end


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
    
    % Put numbers in long format for printing
    format long;
    
    % Read params
    params 

    % random seed
    rng(80085);
    
    % Options for genetic algorithm
    options = gaoptimset('Display','iter','PopulationSize',24,...
    'MutationFcn',@mutationadaptfeasible,...
    'FitnessScalingFcn',@fitscalingrank,'InitialPopulation',pop,'UseParallel','always',...
    'PlotFcns',@gaplotbestf,'EliteCount',0);%,'HybridFcn',{@fmincon,fminconoptions});

    % Network parameter restrictions
    net_lb =  0;
    net_ub =  0.4;
    if rst_type == 'nnt'
        net_lb = 0;
        net_ub = 0;
    end
    
    % Call estimation routine
    [X,fval,exitflag,output,population,scores] = ga(@(X) distance_noprod(X, 0, 1),13, [],[],[],[],[   0.005;  0.01;    6.5;    0.1;  .005; 0.1; 0.1;  0.005; 0.005; 0.5; net_lb; 50; .01], [.5;  1;    10;     1;  0.5;    3;  10; 2; 1; 15; net_ub; 400; 2],[],options);
    
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
end
