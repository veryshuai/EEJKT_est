function simulate(X)
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
    
    % Put numbers in long format for printing
    format long;
    
    % random seed
    rng(80085);
    
    [D,W,error,simulated_data] = distance_noprod(X);
    
    % End timer
    toc
    
    % Close parallel
    matlabpool close
    
    % Save results
    save sim_results;
    
    % Generate Profit Variance Graph
    prof_var(simulated_data);

end
