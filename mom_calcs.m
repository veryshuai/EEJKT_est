% This script contains moment (statistic) calculations for search and learning model

% This script implements the file "system_of_moments_dj2.pdf"

%convert things from cells to matrices
cli_no_mat = cell2mat(cli_no');
sale_h_mat = cell2mat(sale_h');
sale_f_spc = cell2mat(sale_f);
sale_f_mat = cell2mat(sale_f');
ship_f_mat = cell2mat(ship_f');
ship_f_spc = cell2mat(ship_f);
sale_f_mat_count = (sale_f_mat>0);
sh_ann_f_mat = cell2mat(sh_ann_f');
sh_first_yr_dum_mat = cell2mat(sh_first_yr_dum');

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
simulated_data{10} = [cost_f,cost_h];
simulated_data{11} = [prods,succ_prob];
simulated_data{12} = [st_cont,st_ind_cont,cost_vec];

% Create regressors

%% Match based regressors 

% 1. number of clients
cli_no_f = cli_no_mat(:,2:2:end);

% 2 (this one is missing from the initial list) 1 if match exits market
death = zeros(size(sh_ann_f_mat(:)));
death = death.*sh_ann_f_mat(:); % NaNs in sh_ann_f_mat if firm is not exporter
check4nan = zeros(size(death));
check4nan(1:end-1) = isnan(death(2:end)); % is the next element a NaN?
death = death + check4nan; % Only 1 if the next period contains NaN, i.e. a firm death

% 3. 1 if match first year
fy = zeros(size(sh_ann_f_mat(:)));
fy = fy.*sh_ann_f_mat(:); % NaNs in sh_ann_f_mat if firm is not exporter
check4nan = zeros(size(death));
check4nan(2:end) = isnan(death(1:end-1)); % is the last element a NaN?
fy = fy + check4nan; % Only 1 if the last period contains NaN, i.e. match first year 

% 4. Age of match
age = zeros(size(sh_ann_f_mat));
age = age.*sh_ann_f_mat; % NaNs in sh_ann_f_mat if firm is not exporter
%count durations and total exports
for j = 1:size(sh_ann_f_mat,2)
    [first_exp,~] = find(isnan(sh_ann_f_mat(:,j)) == 0,1,'first'); %find the first period of exports
    if isempty(first_exp) == 0 %we want matches that have positive exports 
        for k = 0:size(sh_ann_f_mat,1)-1
            if first_exp + k <= size(sh_ann_f_mat,1) %check to make sure we aren't going past the end of the matrix
                age(first_exp + k,j) = k; % assign age
            end
        end
    end
end



%  
cli_no_f_first = cli_no_mat(1:end-1,:);
cli_no_f_next = cli_no_mat(2:end,:);

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
