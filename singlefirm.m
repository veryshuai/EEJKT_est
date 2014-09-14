function [sh_val_h,sh_val_f,ind,deathmat,ds,sh,act,breakflag,max_mat_violation,match_number_violation,prod_init, cost_vec] = singlefirm(j,max_mat_violation,match_violation,violation,match_number_violation,no_more_rands,breakflag,x_size,Phi_size,z_size,TT,cum_erg_pz,cum_erg_pp,cum_sp_p,th_ind,sp_p,lambda_f_orig,lambda_h_orig,lambda_f_new,lambda_h_new,c_val_h_orig,c_val_f_orig,c_val_h_new,c_val_f_new,burn,delta,d,S,n_size,net_size,Z,Phi,X_f,X_h,actual_h,actual_f,L_b,L_z,L_f,L_h,erg_pz,erg_pp,maxc,max_client_prod,mult_match_max,mms,scale_h,scale_f,de,agg_shocks,cost,F,succ_params)
% This function generates the state matrices for a single firm, this is the main simulation function. 

    % repeatedly reseed in the sequel:
    seed = randi(1e7,6,1);
    
    % Create temp variables (most likely no longer needed)
    exp_inv_temp = 0;
    scalar_temp = 0;
    
    % display(j);
    obin = 20; %current observation index
    ind      = zeros(mms,12); % Event matrix
    succ_prob= zeros(mms,2); % holds success probabilities foreign and home
    deathmat = zeros(mms,1); % Dummy for death of a firm
    ds       = zeros(mms,maxc*2+2); % Demand shocks
    sh       = zeros(mms,maxc*2+2); % Shipment times
    act      = zeros(mms,2); % Dummy for if firm is currently in relationship
    sh_val_h = zeros(mms,3); % Shipment value home (client hotel)
    sh_val_f = zeros(mms,3); % Shipment value foreign (client hotel)
    cost_vec = ones(mms,4) * -1; % flow and fixed, 0 is meaningful
    cost_vec(:,4) = 0; %home fixed cost
    cost_vec(:,2) = 0; %foreign fixed cost
    
    if breakflag == 0
    
        %preallocate
        ind(:,:) = repmat([0,0,0,0,-1,-1,-1,-1,zeros(1,2),-1,-1],mms,1); % holds state
        deathmat(:) = zeros(mms,1); % holds dummy for a firm death
        
        %get Phi's (self productivities)
        get_productivity;
    
        %put in burn period macro shocks (backwards!)
        get_macro_shocks
        
        %resort and fill in
        sort1;
        
        rng(seed(1)); %reseed

        % get foreign matches
        foreign_matches; 
        
        rng(seed(2)); %reseed
        
        % get home matches
        home_matches;
            
        % Small resort 
        [ind(1:obin,:),I] = sortrows((ind(1:obin,:)),1);
        succ_prob(1:obin,:) = succ_prob(I,:);
        deathmat(1:obin,:) = deathmat(I,:);
        cost_vec(1:obin,:) = cost_vec(I,:);
        
        % Read in demand shocks
        demand_shocks;
            
        %now resort and fill in
        sort2;
        
        % fill in the client hotel
        demand_fill_in;
    
        %% Add Shipments
        shipments;
                
        %sort and fill in
        sort3;    
    
        %% Endogenous separation
        end_sep;
        
        %% Active dummy
        active;
    
    end
end
