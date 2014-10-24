function dat = cfs_loop(task,savename,sim_no)
%This function loops simulations of a particular counterfactual, and returns the generated data

    
    dat = cell(sim_no,2);
    
    for k=1:sim_no
        switch task 
            case 1
              simulate(X, 'results/sim_results', 1);
              load('results/cf_sim_results');              
              dat{k,1} = sale_f_mat; dat{k,2} = cell2mat(cli_no'); 
            case 2
              cost_dec_trans(1);
            case 3
              simulate(X, 'results/mac_bump_results', 3)
              load('results/cf_sim_results');              
              dat{k,1} = sale_f_mat; dat{k,2} = cell2mat(cli_no'); 
            case 4
              red_var_trans(1);
            case 5
              search_dec_trans(1);
        end
    end

    save(savename);

end
