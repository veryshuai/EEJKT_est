% This script contains moment (statistic) calculations for search and learning model

% This script implements the file "system_of_moments_dj2.pdf"

%%%%%%%%%%%%%%%%%% Create regressors 

% 1. number of clients
cli_no_mat = cell2mat(cli_no');
cli_no_f = cli_no_mat(:,2:2:end);

% 3. share of first year matches 
match_fy = cell(size(sh_first_yr_dum));
share_fy = cell(size(sh_first_yr_dum));
for j=1:size(sh_first_yr_dum,1) % sum up all first years 
    match_fy{j} = nansum(sh_first_yr_dum{j},2); 
    share_fy{j} = match_fy{j}(1:end-1,:) ./ cli_no{j}(:,2); %extra end row in sh_first_yr_dum was added in st_disc.m for "stacking"
end

% 6. Age of exporter, first and last export year dummies
exp_age = sale_f; exp_first = sale_f; exp_last = sale_f;
%count durations and total export
for l = 1:size(sale_f,1)
    [first_exp,~] = find(sale_f{l} > 0,1,'first'); %find the first period of exports
    [last_exp,~] = find(sale_f{l} > 0,1,'last'); %find the last period of exports
    exp_age{l} = sale_f{l} * NaN; exp_first{l} = sale_f{l} * NaN; exp_last{l} = sale_f{l} * NaN; %this makes the default value NaN rather than zero
    % This loop writes in relationship age
    for k = 0:size(sale_f{1},1)-1
        if first_exp + k <= size(sale_f{1},1) & first_exp + k <= last_exp %check to make sure we aren't going past the end of the matrix
            exp_age{l}(first_exp + k,1) = k + 1; % assign age
            exp_first{l}(first_exp + k,1) = 0; % write in zeros
            exp_last{l}(first_exp + k,1) = 0; % write in zeros
        end
    end
    exp_first{l}(first_exp) = 1; %dummy for first export
    exp_last{l}(last_exp) = 1; %dummy for last export
    if first_exp == 1
        exp_age{l} = sale_f{l} * NaN; % we don't know age if first export in first period -- throw out
        exp_first{l}(first_exp) = NaN; % same logic, don't know first export if observed exporting in first period
    end
    if last_exp == size(sale_f{1},1)
        exp_last{l}(last_exp) = NaN; % we don't know the last export year if observed in last period
    end
end

% 4. Age of match, first and last match period dummies, as well as exporter age
match_age = sh_ann_f; match_first = sh_ann_f; match_last = sh_ann_f; match_exp_age = sh_ann_f;
%count durations and total export
for l = 1:size(sh_ann_f,1)
    for j = 1:size(sh_ann_f{l},2)
        [first_exp,~] = find(isnan(sh_ann_f{l}(:,j)) == 0,1,'first'); %find the first period of exports
        [last_exp,~] = find(isnan(sh_ann_f{l}(:,j)) == 0,1,'last'); %find the last period of exports
        if isempty(first_exp) == 0 %we want matches that have positive exports 
            for k = 0:size(sale_f{1},1)-1
                if first_exp + k <= size(sale_f{1},1) & first_exp + k <= last_exp %check to make sure we aren't going past the end of the matrix
                    match_age{l}(first_exp + k,j) = k + 1; % assign age
                    match_first{l}(first_exp + k,j) = 0; % write in zeros
                    match_last{l}(first_exp + k,j) = 0; % write in zeros
                end
                match_first{l}(first_exp,j) = 1; %dummy for first export
                match_last{l}(last_exp,j) = 1; %dummy for last export
                if first_exp == 1
                    match_age{l}(:,j) = match_age{l}(:,j) * NaN; % we don't know age if first export in first period -- throw out
                    match_first{l}(first_exp,j) = NaN; % same logic, don't know first export if observed exporting in first period
                end
                if last_exp == size(match_age{l},1)
                    match_last{l}(last_exp,j) = NaN; % we don't know the last export year if observed in last period
                end
                match_exp_age{l}(:,j) = [exp_age{l}; NaN]; % add exporter age in the same shape as match data
            end
        end
    end
end

% 5. Average match age by firm
tot_match_age = cell(size(match_age));
avg_match_age = cell(size(match_age));
for j=1:size(match_age,1) % sum up all first years
    tot_match_age{j} = nansum(match_age{j},2);
    avg_match_age{j} = tot_match_age{j}(1:end-1,:) ./ cli_no{j}(:,2); %extra end row in sh_first_yr_dum was added in st_disc.m for "stacking"
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% convert things from cells to matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sale_h_mat = cell2mat(sale_h');
sale_f_spc = cell2mat(sale_f);
sale_f_mat = cell2mat(sale_f');
ship_f_mat = cell2mat(ship_f');
ship_f_spc = cell2mat(ship_f);
sale_f_mat_count = (sale_f_mat>0);
sh_ann_f_mat = cell2mat(sh_ann_f');
sh_first_yr_dum_mat = cell2mat(sh_first_yr_dum');
match_age_mat = cell2mat(match_age'); 
match_first_mat = cell2mat(match_first');
match_last_mat = cell2mat(match_last');
match_exp_age_mat = cell2mat(match_exp_age');

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

%%%%%%%%%% Client Transition Regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Offset client number by one, which will make my life easier during the autoregression
cli_no_next_period = cli_no;
cli_no_for_only = cli_no;
for k = 1:size(cli_no,1)
    cli_no_next_period{k} = [cli_no{k}(2:end,2);NaN];
    cli_no_for_only{k} = cli_no{k}(:,2);
end
cli_coefs = regress(cell2mat(cli_no_next_period),[ones(size(cell2mat(cli_no_next_period))),cell2mat(share_fy),log(cell2mat(cli_no_for_only)),log(cell2mat(cli_no_for_only)).^2,log(cell2mat(avg_match_age)),log(cell2mat(exp_age))]);

exp_death_coefs = regress(cell2mat(exp_last),[ones(size(cell2mat(cli_no_next_period))),cell2mat(share_fy),log(cell2mat(cli_no_for_only)),log(cell2mat(cli_no_for_only)).^2,log(cell2mat(avg_match_age)),log(cell2mat(exp_age))]);

match_death_coefs = regress(match_last_mat(:),[ones(size(match_last_mat(:))),match_first_mat(:),log(sh_ann_f_mat(:)),log(match_age_mat(:)),log(match_exp_age_mat(:))]);

exp_sales_coefs = regress(log(cell2mat(sale_f) - log(cell2mat(cli_no_for_only))),[ones(size(cell2mat(sale_f))),cell2mat(share_fy),cell2mat(cli_no_for_only),log(cell2mat(exp_age))]);

match_ar1_coefs = regress(log(sh_ann_f_mat(2:end))',[ones(size(sh_ann_f_mat(2:end)))',log(sh_ann_f_mat(1:end-1))',match_first_mat(1:end-1)',log(match_exp_age_mat(1:end-1))']);

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
    [dom_ar1_coefs,~,r] = regress(log(lag(:,2)),[ones(size(lag(:,2))),log(lag(:,1))]);
    rMSE = sqrt(sum(r.^2)/size(lag,1));
    mlagdreg = [dom_ar1_coefs(2);rMSE];
catch err
    getReport(err, 'extended') %report error
    mlagdreg = [100;100];
end 

%regression of log(exports) on log(domestic sales) given positive exports
%and positive domestic sales
sale_size = 0;
ind = find(sale_f_mat>sale_size & sale_h_mat>sale_size);
try %allow errors due to no observations etc to be caught
    [exp_dom_coefs,~,r] = regress(log(sale_f_mat(ind)),[ones(size(ind,1),1),log(sale_h_mat(ind))]);
    rMSE = sqrt(sum(r.^2)/size(ind,1));
    mexreg = [exp_dom_coefs(2);rMSE];
catch err
    getReport(err, 'extended') %report error
    mexreg = [100;100];
end

%foreign log log inverse-cdf slope and MSE
cli_no_for_only_mat = cell2mat(cli_no_for_only);
ub = max(cli_no_for_only_mat); %upper bound on client number
inv_cdf = zeros(ub,1);
try %allow errors due to no observations etc to be caught
    for k = 1:ub
        inv_cdf(k) = sum(cli_no_for_only_mat>=k);
    end
    inv_cdf = inv_cdf/sum(inv_cdf);
    [loglog_coefs,~,r] = regress(log(inv_cdf),[ones(ub,1),log((1:ub)'),log((1:ub)').^2]);
    % display(b);
    rMSE = sqrt(sum(r.^2)/ub);
    clidist = [loglog_coefs(2);loglog_coefs(3);rMSE]; %load into moment vector
catch err
    getReport(err, 'extended') %report error
    clidist = [100;100];
end

%average shipments per client per year
cli_no_long = cell2mat(cli_no);
mavship = mean(log(ship_f_spc(ship_f_spc>0)./cli_no_long(ship_f_spc>0,2)));

% Count the number of exporters 
pbexp = sum(sum(sale_f_mat)>1);

