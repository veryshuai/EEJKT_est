function [simulated_data,mavship,loglog_coefs,exp_dom_coefs,dom_ar1_coefs,cli_coefs,exp_death_coefs,match_death_coefs,exp_sales_coefs,match_ar1_coefs] = moms_nocell(mm,lambda_f_orig,lambda_h_orig,lambda_f_new,lambda_h_new,c_val_h_orig,c_val_f_orig,c_val_h_new,c_val_f_new,cf_num,increase,seed)
%This function simulates the model and calculates the moments needed for the distance metric

    % initialize simulated data holder
    simulated_data = cell(11,1);

    % read in parameters
    moms_params
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get Spells for own productivity jumps 
    sp_p = zeros(S,200); %200 possible jumps
    g = 1;
    while (min(sum(sp_p,2))<TT) %loop through jump numbers
        for k = 1:S
          sp_p(k,g) = expinv(rand,1/L_p); 
        end
        g = g + 1;
    end 
    
    %eliminating policy function cell arrrays in favor of multi-dimensional matrices gives the simulation a drastic speed boost
    [lambda_f_orig, lambda_h_orig, c_val_f_orig, c_val_h_orig]  = moms_decell(lambda_f_orig, lambda_h_orig, c_val_f_orig, c_val_h_orig);
    [lambda_f_new, lambda_h_new, c_val_f_new, c_val_h_new]  = moms_decell(lambda_f_new, lambda_h_new, c_val_f_new, c_val_h_new);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Model Simulation

    %% Get vector of state and time changes
    [break_flag,t] = st_traj_nocell(sp_p,lambda_f_orig,lambda_h_orig,lambda_f_new,lambda_h_new,c_val_h_orig,c_val_f_orig,c_val_h_new,c_val_f_new,burn,delta,d,S,n_size,net_size,Z,Phi,X_f,X_h,actual_h,actual_f,L_b,L_z,L_f,L_h,erg_pz,erg_pp,maxc,max_client_prod,mult_match_max,mms,scale_f,scale_h,eta,TT,cost,F,cf_num,succ_params,seed);
    
    % check for errors in simulation routine
    if break_flag == 0
    
        %% Separate dead firms 
        S_old = S;
        sdead(S_old,t);
    
        % check for errors in separation of dead firms 
        if breakflag == 0
    
            %% Calculate Sales
            [sale_h_cont,sale_f_cont] = sales(scale_f,scale_h,eta,st_ind_cont,S,ds,sh,maxc,Z,Phi,X_h,X_f,cf_num,increase,TT);
    
            %% Discretize state vector into years
            [cli_no,sale_h,sale_f,ship_f,sh_ann_f,sh_first_yr_dum,cost_h,cost_f,succ_prob,prods] = st_disc(st_ind_cont,sale_h_cont,sale_f_cont,S,TT,burn,sh,maxc,sh_val_h,sh_val_f,cost_vec,Phi,succ_prob_vec);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Moments calculations
            try 
                mom_calcs;
                % Reject if number of post-burn in exporters less than 500
                display(['Number of post-burn, ever active exporters is ', num2str(pbexp), '.']);
                if pbexp < 500
                    sim_err %fill in parameters to make solver continue
                end
            catch 
                display('Error in moment calculation');
                sim_err %fill in parameters to make solver continue
            end

            % save results 
            if cf_num > 0 & cf_num < 6
                save('results/cf_sim_results') %for plotting counterfactuals
            elseif cf_num == 6
                save('results/val_sim_results') %for calculating the value of the network
            elseif cf_num == 7
                save('results/no_learning_sim_results') %for calculating no learning sales
            elseif cf_num == 8
                save('results/boot_firm_dat') %for bootstrap standard error calculation
            end
    
        else
            %simulation error
            sim_err %fill in parameters to make solver continue
        end
    else
        %simulation error
        sim_err %fill in parameters to make solver continue
    end

    %free memory
    clearvars -except simulated_data mavship loglog_coefs exp_dom_coefs dom_ar1_coefs cli_coefs exp_death_coefs match_death_coefs exp_sales_coefs match_ar1_coefs

end %end function
