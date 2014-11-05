function [D,W,error,simulated_data] = distance_noprod(X, cf_num, seed)
% This function takes a set of parameters and returns the distance between
% the moments in the data and the moments generated by the parameter set.


    % Long format for printing
    format long;

    % Reset random number generator
    if seed
        rng(80085);
    end
    
    % Start timer
    tic
    
    % Change X vector into parameter names
    X2params;

    % Get policy functions
    make_policy;

    % Reset random number generator
    if seed
        rng(80085);
    end

    % Simulate model and create loss function statistics
    [vtran,hazrate,clidist,mstat,mnumex,mavex,mavship,mreg,mexreg,mexshr,mlagereg,mlagdreg,mdeathreg,simulated_data] = moms_nocell(mm,lambda_f_orig,lambda_h_orig,lambda_f_new,lambda_h_new,c_val_h_orig,c_val_f_orig,c_val_h_new,c_val_f_new,cf_num,increase,seed);
    
    %% Targets
    [Data, W] = target_stats();

    %% Simulated data 
    Model = cat(1,vtran,hazrate,clidist,mnumex,mavship,mavex,mreg(1),mreg(3),mexreg,mexshr,mlagereg,mlagdreg,mdeathreg);
    
    % Construct loss
    try
    
        error = Data-Model;
        D = error'*W*error;
    
        % Check for NaNs
        nanflag = isnan(D); 
        if nanflag>0;
            D = D * 10; 
        end
    
        % Punish for value function non-convergence
        punishment;
        D = D*(1+punishment);
        
        % Print Diagnostics
        mmm = cat(2,Data,Model) %data/model comparison
        X %print out current parameter guess
        D %loss
    
        %Simple unweighted loss
        Old_D = norm(error)/norm(Data)
    
        %Print loop time in minutes
        last_loop_run_time = toc/60
    
    catch err

        % report error
        getReport(err, 'extended')
    
        % If broken for any reason, return high loss
        D = 1e12;
    
    end %end try/catch

    % Free up memory
    clearvars -except D W error simulated_data
    % java.lang.Runtime.getRuntime.gc; %java garbage collector

end %end function
