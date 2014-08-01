% This script is the initial call into the search and learning project
% It gives the user the option of which task he would like to perform

% Query user for desired task
acceptable = 'FALSE'; %Is user input interpretable?
while acceptable == 'FALSE'
    task = input(['Task options are currently' char(10) 'est (estimate the model)' char(10) 'sim (simulate the model once)' char(10) 'std (standard error calculation)' char(10) 'cfs (counterfactuals)' char(10) 'val (calc. value of network)' char(10) 'Which task shall I perform?: '], 's');

    display([char(10) 'Your input was ' num2str(task) '.' char(10)]);

    %Check for validity
    if (task == 'est') | (task == 'sim') | (task == 'std') | (task == 'cfs') | (task == 'val')
        acceptable = 'TRUE ';
    else
        display(['Sorry, I do not understand.  Try again.' char(10)]);
    end
end

% Load parameters
params

% Specify Behavior for each task
if task == 'est'
    calibration_noprod(pop, {});
end
if task == 'sim'
    simulate(X, {});
end
if task == 'std'
    bootstrap(X)
    stderr
end
if task == 'cfs'
    call_cfs
end
if task == 'val'
    simulate(X, 'results/val_results', 6)
end




