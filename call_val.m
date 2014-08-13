% calls value of network calucation

    % access demand elasticity and discount from SetParams_noprod
    % the stuff below just makes SetParams run without error
    params; %get X vector
    X2params; %convert X to parameters
    case_str = 'non'; %set counterfactual number
    SetParams_noprod;

    mybeta = 1 - mm.r; %discount
    de = mm.eta; %demand elasticity

    val_mat = cell(50,1); 
    for k=1:size(val_mat,1)
        simulate(X, 'results/val_results', 6);
        val_mat{k} = calc_val(mybeta,de); %get the value of the network
        save('results/network_val_calc_results')
    end
    sum_val
    % Summarize
     
