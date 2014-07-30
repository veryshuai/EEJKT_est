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

% matlabpool open 3
switch task 
    case 1
      simulate(X, 'sim_results', 1);
      cf_decomposition('orig_decomp')
    case 2
      cost_dec_trans(1);
    case 3
      simulate(X, 'mac_bump_results', 3)
      cf_decomposition('mac_bump_decomp')
      makeplots('mac_bump_decomp', 'mac_bump_subplots.eps')
    case 4
      red_var_trans(1);
    case 5
      search_dec_trans(1);
end
%matlabpool close
