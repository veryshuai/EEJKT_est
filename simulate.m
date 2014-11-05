function simulate(X, varargin)
% This script returns simulated data
% for use in creating plots and counterfactuals

    % Start timer 
    tic

    % only allow a two inputs
    numvarargs = size(varargin, 2);
    if numvarargs > 4
        error('distance_noprod:TooManyInputs', ...
            'allow at most 4 optional inputs');
    end

    % set defaults for optional inputs
    optargs = {'sim_results', 0, 0, 1};

    % overwrite defaults with user input
    optargs(1:numvarargs) = varargin;

    % memorable variable names
    savename = optargs{1}; % name underwhich to save results 
    cf_num = optargs{2}; % which counterfactual (0 = none, 6 = calc value, see call_cfs.m for other definitions)?
    debug = optargs{3}; % debug mode -- turn parallel off
    seed = optargs{4}; % seed random number generator?

    % % Parallel setup
    % if debug == 0
    %     if matlabpool('size')~=2 %if pool not equal to 12, open 12
    %        if matlabpool('size')>0
    %          matlabpool close
    %        end
    %        matlabpool open 2
    %     end
    % end
    
    % Put numbers in long format for printing
    format long;
    
    % random seed
    if seed == 1;
        rng(80085);
    end
    
    if debug == 0
        %parfor k=1 %behavior or randoms is different when in parfor loop
        for k=1 %behavior or randoms is different when in parfor loop
            [D,W,err,simulated_data] = distance_noprod(X, cf_num, seed);
        end
    else
        for k = 1 %debug with no parallel
            [D,W,err,simulated_data] = distance_noprod(X, cf_num, seed);
        end
    end
    
    % End timer
    toc
    
    % Save results
    save(savename);
    
    % Close parallel 
    if matlabpool('size')>0
      matlabpool close
    end
    
    % Generate Profit Variance Graph
    %prof_var(simulated_data);

end
