function simulate(X, varargin)
% This script returns simulated data
% for use in creating plots and counterfactuals

    % only allow a two inputs
    numvarargs = length(varargin);
    if numvarargs > 2
        error('distance_noprod:TooManyInputs', ...
            'allow at most 2 optional inputs');
    end

    % set defaults for optional inputs
    optargs = {'sim_results' 1};

    % overwrite defaults with user input
    optargs(1:numvarargs) = varargin;

    % memorable variable names
    savename = optargs{1} % name underwhich to save results 
    cf_num = optargs{2} % which counterfactual (1 = none, see call.m for other definitions)?

    % Parallel setup
    clc
    if matlabpool('size')~=12 %if pool not equal to 12, open 12
       if matlabpool('size')>0
         matlabpool close
       end
       matlabpool open 12
    end
    
    % Put numbers in long format for printing
    format long;
    
    % random seed
    rng(80085);
    
    [D,W,err,simulated_data] = distance_noprod(X, cf_num);
    
    % End timer
    toc
    
    % Close parallel
    matlabpool close
    
    % Save results
    save(savename);
    
    % Generate Profit Variance Graph
    prof_var(simulated_data);

end
