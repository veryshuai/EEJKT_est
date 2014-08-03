function [vtran,hazrate,clidist,mstat,mnumex,mavex,mavship,mreg,mexreg,mexshr,mlagereg,mlagdreg,mdeathreg,simulated_data] = moms_nocell(mm,lambda_f_orig,lambda_h_orig,lambda_f_new,lambda_h_new,c_val_h_orig,c_val_f_orig,c_val_h_new,c_val_f_new,cf_num,increase)
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate acceptance rates for home and foreign market

    %Initialize the theta distribution parameter vectors
    th0 = zeros(S,1);
    th1 = zeros(S,1);
    th2 = zeros(S,1);

    %Get random draws of theta for each firm
    for k = 1:S
        th0(k) = betainv(rand,ag,bg);
        th1(k) = betainv(rand,ah,bh);
        th2(k) = betainv(rand,af,bf);
    end
    th_draw   = cat(2,th0,th1,th2);
    
    % Read in the theta grid values 
    th    = cell(3,1);
    th{1} = theta0;
    th{2} = theta1;
    th{3} = theta2;
    
    % Find the closest grid values to each firms random draw
    indx1 = zeros(S,3);
    theta = zeros(S,3);
    for j = 1:3  % Map theta0, theta1, and theta2 draws onto grid.  
      for s = 1:S
        dif            = abs(th_draw(s,j)-th{j}); 
        [~,indx1(s,j)] = min(dif);
        theta(s,j)     = th{j}(indx1(s,j));
     end
    end
    
    % Get common vs independent components in success probabilities
    mu_h = myalpha*theta(:,1)+(1-myalpha)*theta(:,2); %true home success probability
    mu_f = myalpha*theta(:,1)+(1-myalpha)*theta(:,3); %true foreign success probability
    
    %eliminating policy function cell arrrays in favor of multi-dimensional matrices gives the simulation a drastic speed boost
    [lambda_f_orig, lambda_h_orig, c_val_f_orig, c_val_h_orig]  = moms_decell(lambda_f_orig, lambda_h_orig, c_val_f_orig, c_val_h_orig);
    [lambda_f_new, lambda_h_new, c_val_f_new, c_val_h_new]  = moms_decell(lambda_f_new, lambda_h_new, c_val_f_new, c_val_h_new);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Model Simulation

    %% Get vector of state and time changes
    [st_ind_cont,st_cont,ds,sh,act,break_flag,deathmat,sh_val_h,sh_val_f,cprod] = st_traj_nocell(indx1,mu_h,mu_f,sp_p,lambda_f_orig,lambda_h_orig,lambda_f_new,lambda_h_new,c_val_h_orig,c_val_f_orig,c_val_h_new,c_val_f_new,burn,delta,d,S,n_size,net_size,Z,Phi,X_f,X_h,actual_h,actual_f,L_b,L_z,L_f,L_h,erg_pz,erg_pp,maxc,max_client_prod,mult_match_max,mms,scale_f,scale_h,eta,TT);
    
    % check for errors in simulation routine
    if break_flag == 0
    
        %% Separate dead firms 
        S_old = S;
        [st_cont,st_ind_cont,S,ds,sh,breakflag,sh_val_h,sh_val_f] = sdead(st_cont,st_ind_cont,S,ds,sh,deathmat,sh_val_h,sh_val_f);
    
        % check for errors in separation of dead firms 
        if breakflag == 0
    
            %% Calculate Sales
            [sale_h_cont,sale_f_cont] = sales(scale_f,scale_h,eta,st_ind_cont,S,ds,sh,maxc,Z,Phi,X_h,X_f,cf_num,increase,TT);
    
            %% Discretize state vector into years
            [cli_no,sale_h,sale_f,ship_f,sh_ann_f,sh_first_yr_dum] = st_disc(st_ind_cont,sale_h_cont,sale_f_cont,S,TT,burn,sh,maxc,sh_val_h,sh_val_f);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Moments calculations
    
            %convert things from cells to matrices
            cli_no_mat = cell2mat(cli_no);
            sale_h_mat = cell2mat(sale_h');
            sale_f_spc = cell2mat(sale_f);
            sale_f_mat = cell2mat(sale_f');
            ship_f_mat = cell2mat(ship_f');
            ship_f_spc = cell2mat(ship_f);
            sale_f_mat_count = (sale_f_mat>0);
            sh_ann_f_mat = cell2mat(sh_ann_f');
            sh_first_yr_dum_mat = cell2mat(sh_first_yr_dum');
            prods = cell2mat(cprod);

            %pack simulation info up to return it
            simulated_data{1} = cli_no_mat;
            simulated_data{2} = sale_h_mat;
            simulated_data{3} = sale_f_spc;
            simulated_data{4} = sale_f_mat;
            simulated_data{5} = ship_f_mat;
            simulated_data{6} = ship_f_spc;
            simulated_data{7} = sale_f_mat_count;
            simulated_data{8} = sh_ann_f_mat;
            simulated_data{9} = sh_first_yr_dum_mat;
            simulated_data{10} = [prods,mu_f,mu_h];
            simulated_data{11} = [st_cont,st_ind_cont];
    
            %client transition counts
            trans_f = zeros(4);
            for j = 1:S
              trans_f_inter = sparse(min(cli_no{j}(1:end-1,2),3)+1,min(cli_no{j}(2:end,2),3)+1,1,4,4); %interim count
              trans_f = trans_f + trans_f_inter;
            end 
            trans_f = full(trans_f)./repmat(sum(trans_f,2),1,size(trans_f,2));
            vtran = [trans_f(2,1:4)';trans_f(3,1:4)'];  
            %vtran is the transitions from (one,two) to (zero,one,two,three or more) foreign clients
    
            %count export durations and total exports
            durcount = zeros(7,1); %number of exporters of age k 
            totex = zeros(7,1); %total export sales in age k of being an exporter
            l_tot = zeros(7,1); %sum of log export sales in age k of being an exporter
            for j = 1:S
                [first_exp,~] = find(sale_f_mat_count(:,j),1,'first'); %find the first period of exports
                if isempty(first_exp) == 0 && first_exp ~= 1 %we want firms that have positive exports, and we don't count the first cohort (because they might have already been exporting.
                    for k = 1:7
                        if first_exp + k-1 <= size(sale_f_mat_count,1) %check to make sure we aren't going past the end of the matrix
                            if sale_f_mat(first_exp+k-1,j)>0 %don't count inactive exporters
                                durcount(k) = durcount(k) + 1;
                                totex(k) = totex(k) + sale_f_mat(first_exp+k-1,j); 
                                l_tot(k) = l_tot(k) + log(sale_f_mat(first_exp+k-1,j));
                            end
                        end
                    end
                end
            end
            
            %log average sales in year k of being an exporter
            mavex = l_tot(1:5)./durcount(1:5);
            
            %export separation hazard, conditional on survival
            mnumex = 1-durcount(2:6)./durcount(1:5); 
            mnumex(isnan(mnumex) == 1) = 5; %mnumex can be 0/0, so check for NaN to be safe
            
            %match separation hazard, conditional on survival
            active_to_k = zeros(5,1);
            active_to_gtk = zeros(5,1);
            for j=1:S_old
                for k=1:5
                    active_to_k(k) = active_to_k(k) + sum(act{j}(:,end)==k);
                    active_to_gtk(k) = active_to_gtk(k) + sum(act{j}(:,end)>=k); 
                end
            end
            hazrate = active_to_k./active_to_gtk;
            hazrate(isnan(hazrate) == 1) = 5; %hazrate can be 0/0, so check for NaN to be safe
            
            %foreign log log inverse-cdf slope and MSE
            ub = max(cli_no_mat(:,2)); %upper bound on client number
            inv_cdf = zeros(ub,1);
            if ub > 0;
                for k = 1:ub
                    inv_cdf(k) = sum(cli_no_mat(:,2)>=k);
                end
                inv_cdf = inv_cdf/sum(inv_cdf);
                [b,~,r] = regress(log(inv_cdf),[ones(ub,1),log((1:ub)'),log((1:ub)').^2]);
                % display(b);
                rMSE = sqrt(sum(r.^2)/ub);
                clidist = [b(2);b(3);rMSE];
            else
                clidist = [100;100];
                display('ERROR: Nothing to regress moms.m client distribution');
            end
            
            %average log domestic sales, log foreign sales, and standard deviations.
            sale_size = 0; %in case we want to add size restriction to match establishment survey rule
            mstat = zeros(2,1);
            mstat(1) = mean(log(sale_f_mat(sale_f_mat>sale_size)));
            mstat(2) = sqrt(mean((log(sale_f_mat(sale_f_mat>sale_size))-mstat(1)).^2));
            
            %regression of log(exports) on log(domestic sales) given positive exports
            %and positive domestic sales
            ind = find(sale_f_mat>sale_size & sale_h_mat>sale_size);
            if isempty(ind) == 0
                [b,~,r] = regress(log(sale_f_mat(ind)),[ones(size(ind,1),1),log(sale_h_mat(ind))]);
                rMSE = sqrt(sum(r.^2)/size(ind,1));
                mexreg = [b(2);rMSE];
            else
                mexreg = [100;100];
                display('ERROR: Nothing to regress moms.m exports on domestic sales regression');
            end
            
            %share of exporting plants among active plants
            mexshr = sum(sale_f_mat>0)/sum(sale_f_mat>0 | sale_h_mat>0);
            
            %regression of log sales per client on log number of clients, conditonal on
            %positive sales.
            spc = sale_f_spc./cli_no_mat(:,2);
            ind = find(isnan(spc) ==0 & spc>0);
            if isempty(ind) == 0 && rank([ones(size(ind,1),1),log(cli_no_mat(ind,2)),log(cli_no_mat(ind,2)).^2]) == 3 %make sure we have full rank
                [b,~,r] = regress(log(spc(ind)),[ones(size(ind,1),1),log(cli_no_mat(ind,2)),log(cli_no_mat(ind,2)).^2]);              
                rMSE = sqrt(sum(r.^2)/size(ind,1));
                mreg = [b(2);b(3);rMSE];
            else
                mreg = [100;100;100];
                display('ERROR: Nothing to regress moms.m sales-per-client regression');
            end
            if isnan(mreg) ~= [0;0;0]
                mreg = [100;100;200];
                display('ERROR: nans in mreg');
            end
            
            %average shipments per client per year
            mavship = mean(ship_f_spc(cli_no_mat(:,2)>0)./cli_no_mat(cli_no_mat(:,2)>0,2));
            
            %regression of log domestic sales on log lagged domestic sales
            lag1 = zeros(9*S,2);
            lag2 = zeros(8*S,2);
            for k = 1:9
              lag1((k-1)*S+1:k*S,:) = sale_h_mat(2*(k-1)+1:2*(k-1)+2,:)';
            end
            for k = 1:8
              lag2((k-1)*S+1:k*S,:) = sale_h_mat(2*(k-1)+2:2*(k-1)+3,:)';
            end
            lag = [lag1;lag2];
            lag = lag((lag(:,1)>0 & lag(:,2)>0),:);
            if isempty(lag) == 0
                [b,~,r] = regress(log(lag(:,2)),[ones(size(lag(:,2))),log(lag(:,1))]);
                rMSE = sqrt(sum(r.^2)/size(lag,1));
                mlagdreg = [b(2);rMSE];
            else
               mlagdreg = [100;100];
               display('ERROR: Nothing to regress moms.m domestic sales ar1 regression');
            end 
            
            % Count the number of exporters 
            pbexp = sum(sum(sale_f_mat)>1);
            
            % match level ar1
            try
                [b,~,r] = regress(log(sh_ann_f_mat(2:end))',[ones(size(sh_ann_f_mat(2:end)))',log(sh_ann_f_mat(1:end-1))',sh_first_yr_dum_mat(1:end-1)']);
                rMSE = sqrt(sum(r(isnan(r)==0).^2)/size(sh_ann_f_mat(isnan(sh_ann_f_mat)==0),1));
                mlagereg = [b(2);b(3);rMSE];
            catch err
               mlagereg = [100;100];
               display('ERROR: Nothing to regress moms.m match ar1 regression');
            end 
            
            %death regression
            % first a bit of manipulation -- creat a matrix of annual sales of same
            % size as sales by client matrix
            tot_sales_cell = cell(S,1);
            %ugly for-loop to get annual totals in same size as sh_ann_f_mat
            for j = 1:S
                tot_sales_cell{j} = repmat(sum(max(sh_ann_f{j},0),2),1,size(sh_ann_f{j},2));
            end
            tot_sales_mat = cell2mat(tot_sales_cell');
            
            %now we need a dummy which is one when there is a death, and zero when no
            %death
            death = zeros(size(sh_ann_f_mat(:)));
            death = death.*sh_ann_f_mat(:); % put nans in correct places
            check4nan = zeros(size(death));
            check4nan(1:end-1) = isnan(death(2:end)); % is the next element a NaN?
            death = death + check4nan; % this gives us what we want
            
            %regression
            rhs = [ones(size(death)),sh_first_yr_dum_mat(:),log(tot_sales_mat(:))];
            try
                [b,~,r] = regress(death,rhs);
                rMSE = sqrt(sum(r(isnan(r)==0).^2)/size(death(isnan(death)==0),1));
                mdeathreg = [b;rMSE];
            catch err
                mdeathreg = ones(4,1) * 100;
                display('funny business in death regression')
            end
            
            display(['Number of post-burn, ever active exporters is ', num2str(pbexp), '.']);

            % save results 
            if cf_num > 0 & cf_num < 6
                save('results/cf_sim_results') %for plotting counterfactuals
            elseif cf_num == 6
                save('results/val_sim_results') %for calculating the value of the network
            end
    
            % Reject if number of post-burn in exporters less than 500
            if pbexp < 500
                sim_err %fill in parameters to make solver continue
            end
        else
            %simulation error
            sim_err %fill in parameters to make solver continue
        end
    else
        %simulation error
        sim_err %fill in parameters to make solver continue
    end

end %end function
