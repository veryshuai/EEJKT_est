function bootstrap(X)
%This function takes parameters solved for in the "calibration" optimization loop, generates a bunch of firms, and bootstraps standard errors for the moments which are used in the indirect inference part of the optimization program.

    %number of bootstrap draws
    boot_num = 2000;
    
    if matlabpool('size')~=7
      if matlabpool('size')>0
        matlabpool close
      end
      matlabpool open 7
    end
    
    % load data
    load results/boot_firm_dat
    
    display(['We created a total of ', num2str(S), 'firms.']);
    
    %% Calculate the moments, with resampling
    
    %how many exporters? from Table 1, summing over all non-affiliated cohorts we get a total of 20,990
    %exporters.  Of course, some are only observed for a single year, and the first year picks up all
    %who were already exporting.
    exp_num = 20990;
    S = exp_num;
    
    %clear unnecessary variables
    clearvars -except act cli_no sale_h sale_f ship_f S S_old burn X boot_num exp_num sh_ann_f sh_first_yr_dum
    
    %eliminate inactive firms
    sale_f_mat      = cell2mat(sale_f');
    sale_h_mat      = cell2mat(sale_h');
    ind             = logical(sum(sale_f_mat)==0 & sum(sale_h_mat) == 0)';	
    cli_no          = cli_no(~ind);
    sale_h          = sale_h(~ind);
    sale_f          = sale_f(~ind);
    ship_f          = ship_f(~ind);
    sh_ann_f        = sh_ann_f(~ind);
    sh_first_yr_dum = sh_first_yr_dum(~ind);
    
    %create empty matrix to hold moments
    moments = zeros(41,boot_num);
    
    % Reset random seed
    rng(37738); %belle upside down
    
    % MAIN MOMENT CREATION LOOP (COPIED FROM MOMS.M)
    parfor ctr = 1:boot_num
      S_old = 10000;
      display(ctr);
      
      %draw clients
      ind = randi(size(sale_f,1),exp_num,1);
      boot_cli_no = cli_no(ind);
      boot_sh_ann_f = sh_ann_f(ind);
      
      %convert things from cells to matrices
      cli_no_mat = cell2mat(boot_cli_no);
      sale_h_mat = cell2mat(sale_h(ind)');
      sale_f_spc = cell2mat(sale_f(ind));
      sale_f_mat = cell2mat(sale_f(ind)');
      ship_f_mat = cell2mat(ship_f(ind)');
      ship_f_spc = cell2mat(ship_f(ind));
      sale_f_mat_count = (sale_f_mat>0);
      sh_ann_f_mat = cell2mat(sh_ann_f(ind)');
      sh_first_yr_dum_mat = cell2mat(sh_first_yr_dum(ind)');
      sh_ann_f_death = sh_ann_f(ind); %for the death calculation
    
      %client transition counts
      trans_f = zeros(4);
      for j = 1:S
        trans_f_inter = sparse(min(boot_cli_no{j}(1:end-1,2),3)+1,min(boot_cli_no{j}(2:end,2),3)+1,1,4,4); %interim count
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
      try %allow errors due to no observations etc to be caught
          for k = 1:ub
              inv_cdf(k) = sum(cli_no_mat(:,2)>=k);
          end
          inv_cdf = inv_cdf/sum(inv_cdf);
          [b,~,r] = regress(log(inv_cdf),[ones(ub,1),log((1:ub)'),log((1:ub)').^2]);
          % display(b);
          rMSE = sqrt(sum(r.^2)/ub);
          clidist = [b(2);b(3);rMSE]; %load into moment vector
      catch err
          getReport(err, 'extended') %report error
          clidist = [100;100];
      end
      
      %average log domestic sales, log foreign sales, and standard deviations.
      sale_size = 0; %in case we want to add size restriction to match establishment survey rule
      mstat = zeros(2,1);
      mstat(1) = mean(log(sale_f_mat(sale_f_mat>sale_size)));
      mstat(2) = sqrt(mean((log(sale_f_mat(sale_f_mat>sale_size))-mstat(1)).^2));
      
      %regression of log(exports) on log(domestic sales) given positive exports
      %and positive domestic sales
      ind = find(sale_f_mat>sale_size & sale_h_mat>sale_size);
      try %allow errors due to no observations etc to be caught
          [b,~,r] = regress(log(sale_f_mat(ind)),[ones(size(ind,1),1),log(sale_h_mat(ind))]);
          rMSE = sqrt(sum(r.^2)/size(ind,1));
          mexreg = [b(2);rMSE];
      catch err
          getReport(err, 'extended') %report error
          mexreg = [100;100];
      end
      
      %share of exporting plants among active plants
      mexshr = sum(sale_f_mat>0)/sum(sale_f_mat>0 | sale_h_mat>0);
      
      %regression of log sales per client on log number of clients, conditonal on
      %positive sales.
      spc = sale_f_spc./cli_no_mat(:,2);
      ind = find(isnan(spc) ==0 & spc>0);
      try %allow errors due to no observations etc to be caught
          [b,~,r] = regress(log(spc(ind)),[ones(size(ind,1),1),log(cli_no_mat(ind,2)),log(cli_no_mat(ind,2)).^2]);              
          rMSE = sqrt(sum(r.^2)/size(ind,1));
          mreg = [b(2);b(3);rMSE];
      catch err
          getReport(err, 'extended') %report error
          mreg = [100;100;100];
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
      try
          [b,~,r] = regress(log(lag(:,2)),[ones(size(lag(:,2))),log(lag(:,1))]);
          rMSE = sqrt(sum(r.^2)/size(lag,1));
          mlagdreg = [b(2);rMSE];
      catch err
          getReport(err, 'extended') %report error
          mlagdreg = [100;100];
      end 
      
      % Count the number of exporters 
      pbexp = sum(sum(sale_f_mat)>1);
      
      % match level ar1
      try
          [b,~,r] = regress(log(sh_ann_f_mat(2:end))',[ones(size(sh_ann_f_mat(2:end)))',log(sh_ann_f_mat(1:end-1))',sh_first_yr_dum_mat(1:end-1)']);
          rMSE = sqrt(sum(r(isnan(r)==0).^2)/size(sh_ann_f_mat(isnan(sh_ann_f_mat)==0),1));
          mlagereg = [b(2);b(3);rMSE];
      catch err
          getReport(err, 'extended') %report error
          mlagereg = [100;100];
      end 
      
      %death regression
      % first a bit of manipulation -- creat a matrix of annual sales of same
      % size as sales by client matrix
      tot_sales_cell = cell(S,1);
      %ugly for-loop to get annual totals in same size as sh_ann_f_mat
      for j = 1:S
          tot_sales_cell{j} = repmat(sum(max(boot_sh_ann_f{j},0),2),1,size(boot_sh_ann_f{j},2));
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
      
      moments(:,ctr) = cat(1,vtran,hazrate,clidist,mnumex,mavship,mavex,mreg(1),mreg(3),mexreg,mexshr,mlagereg,mlagdreg,mdeathreg);
    end
    
    bcov = cov(moments');
    
    save results/bcov;

end
