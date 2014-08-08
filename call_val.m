% calls value of network calucation

    beta = 0.95; % annual discount
    val_mat = cell(20,1); 
    for k=1:20
        simulate(X, 'results/val_results', 6);
        val_mat{k} = calc_val(beta); %get the value of the network
        save('results/network_val_calc_results')
    end
    sum_val
    % Summarize
     
