function [cind,cst,cds,csh,cact,breakflag,cdeathmat,csh_val_h,csh_val_f,cprod] = st_traj_nocell(th_ind,mu_h,mu_f,sp_p,lambda_f,lambda_h,c_val_h,c_val_f,burn,delta,d,S,n_size,net_size,Z,Phi,X_f,X_h,actual_h,actual_f,L_b,L_z,L_f,L_h,erg_pz,erg_pp,maxc,max_client_prod,mult_match_max,mms,scale_f,scale_h,de,esT)
    %This function does the simulation needed to calculate moments
    
    breakflag = 0; %this flag goes to one if there is a maximum matrix violation, and allows us to stop the loop  
    
    match_violation = 0; %this counts number of match per period violations
    max_mat_violation = 0; %counts number of matrix size violations
    violation = 0; %will count number of demand shock change violations
    match_number_violation = 0; %counts number of maximum match number violations
    no_more_rands = 0; %counts number of times end of random vector reached 
    
    %% SETTING UP
    
    % Initialize cells
    cind         = cell(S,1);
    cst          = cell(S,1);
    cds          = cell(S,1);
    csh          = cell(S,1);
    cact         = cell(S,1);
    cdeathmat    = cell(S,1);
    csh_val_h    = cell(S,1);
    csh_val_f    = cell(S,1);
    cprod        = cell(S,1);
    
    % cumulative versions of ergodic distributions and waiting times
    cum_erg_pz = cumsum(erg_pz);
    cum_erg_pp = cumsum(erg_pp);
    cum_sp_p = cumsum(sp_p,2);
    
    %find middle state
    Phi_mid = (size(Phi,1)+1)/2;
    X_f_mid = (size(X_f,1)+1)/2;
    Z_mid = (size(Z,1) +1)/2;
    
    x_size = size(X_f,1)-X_f_mid;
    Phi_size = size(Phi,1)-Phi_mid;
    z_size = size(Z,1)-Z_mid;
    
    %get firm specific simulation times
    TT = esT;
    
    %reset random seed
    rng(80085);
    
    %create common aggregate shocks (for burn in)
    agg_shocks = rand(10000,1);
    
    %give each firm a random seed
    seeds = randi(1e7,S,1);
    
    %% EXOGENOUS TRAJECTORIES (MACRO STATES AND SELF PRODUCTIVITY)
    for j = 1:S
        display(j);
        rng(seeds(j));
        nonactives = 0;
        if breakflag == 0
            notempty = 0;
            counter = 0;
            while notempty == 0
                if counter<3
                    [sh_val_h,sh_val_f,ind,deathmat,ds,sh,act,breakflag,max_mat_violation,match_number_violation,prod_init] = singlefirm(j,max_mat_violation,match_violation,violation,match_number_violation,no_more_rands,breakflag,x_size,Phi_size,z_size,TT,cum_erg_pz,cum_erg_pp,cum_sp_p,th_ind,mu_h,mu_f,sp_p,lambda_f,lambda_h,c_val_h,c_val_f,burn,delta,d,S,n_size,net_size,Z,Phi,X_f,X_h,actual_h,actual_f,L_b,L_z,L_f,L_h,erg_pz,erg_pp,maxc,max_client_prod,mult_match_max,mms,scale_h,scale_f,de,agg_shocks);
                    burn_ind = find((ind(:,1))>=burn,1,'first');
                    if isempty(burn_ind)==0
                        notempty = (sum(sum(sh(burn_ind:end,:))));
                    end
                    counter = counter + 1;
                else
                    notempty = 1;
                    nonactives = nonactives + 1;
                    if nonactives > 2000
                        display('ERROR: Too many trials with no observations');
                        breakflag = 1;
                    end
                end
            end
        lt = find((ind(:,1))'>0,1,'last');
        ft = find((ind(:,1))'>0,1,'first')+1;
    
        % read into sparses
        cind{j}       = sparse(ind(ft:lt,:));
        cds{j}        = sparse(ds(ft:lt,:));
        csh{j}        = sparse(sh(ft:lt,:));
        cact{j}       = sparse(act(ft:lt,:));
        cdeathmat{j}  = sparse(deathmat(ft:lt,:));
        csh_val_h{j}  = sparse(sh_val_h(ft:lt,:));
        csh_val_f{j}  = sparse(sh_val_f(ft:lt,:));
        cprod{j}      = sparse(prod_init);
        end
    end
    
    display(['A total of ', num2str(match_number_violation),' maximum match violations']);
    display(['A total of ', num2str(max_mat_violation),' matrix size violations']);
    
    if breakflag == 1
        display('WARNING: Broke out of loop! Results not reliable.')
    end
end
