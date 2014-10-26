% This script calls counterfactuals (need to do this due to parallel issues)
% In particular, to get the same random draws as in estimation we need to be running in parallel, but I get a transparency error if I directly put parallel into the scripts 

% Query user for desired task
acceptable = 'FALSE'; %Is user input interpretable?
while acceptable == 'FALSE'
    task = input([char(10) 'Please enter the number of the counterfactual you would like to run. Currently we have:' char(10) '1. Simple simulation' char(10) '2. Decrease in fixed cost of maintaining a relationship' char(10) '3. Increase in macro shock' char(10) '4. Reduction in variance of macro shock' char(10) '5. Decrease in search cost' char(10)]);

    display([char(10) 'Your input was ' num2str(task) '.' char(10)]);

    %Check for validity
    if (task == 1) | (task == 2) | (task == 3) | (task == 4) | (task == 5)
        acceptable = 'TRUE ';
    else
        display(['Sorry, I do not understand.  Try again.' char(10)]);
    end
end

% Query user for desired number of simulations to average with random macro shocks
sim_no = input([char(10) 'Please enter the desired number of random macro shock simulations you would like to average (default 10):' char(10)]);

display([char(10) 'Your input was ' num2str(sim_no) '.' char(10)]);

%Check user input for validity
if isnumeric(sim_no) == 0 | sim_no ~= round(sim_no) | sim_no < 0
    display(['Sorry, I do not understand.  Using default value of 10']);
    sim_no = 10;
end

% Counterfactual loop (many macro shocks) and plot
% matlabpool open 3
switch task 
    case 1
        cfs_loop(task,'results/orig',sim_no);
        cf_decomposition('results/orig','results/orig_decomp');
    case 2
        pass;
    case 3
        cfs_loop(task,'results/mac_bump',sim_no);
        cf_decomposition('results/mac_bump','results/mac_bump_decomp');
        makeplots('results/mac_bump_decomp', 'results/mac_bump_subplots.eps')
    case 4
        pass;
    case 5
        pass;
end
%matlabpool close
