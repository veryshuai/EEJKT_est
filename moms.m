function [vtran,hazrate,clidist,mstat,mnumex,mavex,mavship,mreg,mexreg,mexshr]...
= moms(mm,c_val_h,c_val_f,lambda_f,lambda_h,rv)
%This function calculates the moments needed for the distance metric

burn    = mm.burn;    %number of burn in periods
esT     = 18 + burn;  % number of ergodic state periods to be simulated
scale_f = mm.scale_f;  %US scale parameter
scale_h = mm.scale_h;  %Colombia scale parameter
eta     = mm.eta;      %demand elasticity
alpha   = mm.alpha;    %success function parameter
ah      = mm.ah;  %home Beta parameter
bh      = mm.bh;  %home Beta parameter  
af      = mm.af;  %foreign Beta parameter
bf      = mm.bf;  %foreign Beta parameter
S       = mm.S;   %number of firms

ag    = mm.ag;     % true common (theta0) beta parameter
bg    = mm.bg;     % true common (theta0) beta parameter
theta0 = mm.theta0;  % vector of global firm effects
theta1 = mm.theta2;  % vector of home market firm effects
theta2 = mm.theta2;  % vector of foreign market firm effects

TT          = esT;           % Number of time periods to be simulated

L_p         = mm.L_p;      %arrival rate for jumps in self productivity
L_h         = mm.L_h;
L_f         = mm.L_f;

d           = mm.d;  %exogenous death rate

maxc        = mm.maxc;      %maximum allowable number of clients

%% Get Spells for Own productivity jumps (need to do this every iteration since we are estimating the jump probabilities)

sp_p = [];
while (isempty(sp_p)==1) || (min(sum(sp_p,2))<TT)
    sp_p = [sp_p,zeros(S,1)];
    for k = 1:S
      sp_p(k,end) = expinv(rv(ri),1/L_p); %note that the parameter is 1/rate due to how matlab calls exponential distribution
      ri = mod(ri,rs)+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate acceptance rates for home and foreign market
th0 = zeros(S,1);
th1 = zeros(S,1);
th2 = zeros(S,1);
for k = 1:S
    th0(k) = betainv(rv(ri),ag,bg);
    ri = mod(ri,rs)+1;
    th1(k) = betainv(rv(ri),ah,bh);
    ri = mod(ri,rs)+1;
    th2(k) = betainv(rv(ri),af,bf);
    ri = mod(ri,rs)+1;
end
th_draw   = cat(2,th0,th1,th2);

th    = cell(3,1);
th{1} = theta0;
th{2} = theta1;
th{3} = theta2;

indx1 = zeros(S,3);
theta = zeros(S,3);

for j = 1:3  % Map theta0, theta1, and theta2 draws onto grid.  
  for s = 1:S
    dif            = abs(th_draw(s,j)-th{j}); 
    [~,indx1(s,j)] = min(dif);
    theta(s,j)     = th{j}(indx1(s,j));
 end
end

mu_h = alpha*theta(:,1)+(1-alpha)*theta(:,2); %true home success probability
mu_f = alpha*theta(:,1)+(1-alpha)*theta(:,3); %true foreign success probability

%% Get vector of state and time changes
[st_ind_cont,st_cont,ds,sh,act,break_flag,deathmat] = st_traj(indx1,mu_h,mu_f,sp_p,lambda_f,lambda_h,c_val_h,c_val_f,mm,ri,rv,rs); %indx1 is a theta index (firm;general,home,foreign)

if break_flag == 0

  clearvars -except st_ind_cont st_cont eta scale_f scale_h S TT burn ds sh act maxc deathmat %clear memory of unecessary stuff

  %% Separate dead firms
  S_old = S;
  [st_cont,st_ind_cont,S,ds,sh,breakflag] = sdead(st_cont,st_ind_cont,S,ds,sh,deathmat);

  if breakflag == 0
    %% Sales
    [sale_h_cont,sale_f_cont] = sales(scale_f,scale_h,eta,st_cont,S,ds,sh,maxc);

    %% Discretize state vector
    [cli_no,sale_h,sale_f,ship_f] = st_disc(st_ind_cont,sale_h_cont,sale_f_cont,S,TT,burn,sh,maxc);

    clearvars -except act cli_no sale_h sale_f ship_f S S_old burn %clear memory of unnecessary stuff

    %% Moments calculations

    %convert things from cells to matrices
    cli_no_mat = cell2mat(cli_no);
    sale_h_mat = cell2mat(sale_h');
    sale_f_spc = cell2mat(sale_f);
    sale_f_mat = cell2mat(sale_f');
    ship_f_mat = cell2mat(ship_f');
    ship_f_spc = cell2mat(ship_f);
    sale_f_mat_count = (sale_f_mat>0);

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
%mavex = l_tot(1:5)./durcount(1:5);
mavex = log(totex(1:5)./durcount(1:5));

%export separation hazard, conditional on survival
mnumex = 1-durcount(2:6)./durcount(1:5); 
mnumex(isnan(mnumex) == 1) = 1; %mnumex can be 0/0, so check for NaN to be safe

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
hazrate(isnan(hazrate) == 1) = 1; %hazrate can be 0/0, so check for NaN to be safe

%foreign log log inverse-cdf slope and MSE
ub = max(cli_no_mat(:,2)); %upper bound on client number
inv_cdf = zeros(ub,1);
if ub > 0;
    for k = 1:ub
        inv_cdf(k) = sum(cli_no_mat(:,2)>=k);
    end
    inv_cdf = inv_cdf/sum(inv_cdf);
    [b,~,r] = regress(log(inv_cdf),[ones(ub,1),log((1:ub)')]);
    rMSE = sqrt(sum(r.^2)/ub);
    clidist = [b(2);rMSE];
else
    clidist = [10;10];
    display('ERROR: Nothing to regress moms.m client distribution');
end

%average log domestic sales, log foreign sales, and standard deviations.
sale_size = 0; %in case we want to add size restriction to match establishment survey rule
mstat = zeros(2,1);
mstat(1) = log(mean(sale_f_mat(sale_f_mat>sale_size)));
mstat(2) = sqrt(exp(mstat(1))^(-2) * mean((sale_f_mat(sale_f_mat>sale_size)-exp(mstat(1))).^2));

%regression of log(exports) on log(domestic sales) given positive exports
%and positive domestic sales
ind = find(sale_f_mat>sale_size & sale_h_mat>sale_size);
if isempty(ind) == 0
    [b,~,r] = regress(log(sale_f_mat(ind)),[ones(size(ind,1),1),log(sale_h_mat(ind))]);
    rMSE = sqrt(sum(r.^2)/size(ind,1));
    mexreg = [b(1);b(2);rMSE];
else
    mexreg = [10;10;10];
    display('ERROR: Nothing to regress moms.m exports on domestic sales regression');
end

%share of exporting plants among active plants
mexshr = sum(sale_f_mat>0)/sum(sale_f_mat>0 | sale_h_mat>0);

%regression of log sales per client on log number of clients, conditonal on
%positive sales.
spc = sale_f_spc./cli_no_mat(:,2);
[ind,~] = find(isnan(spc) ==0 & spc>0);
if isempty(ind) == 0 && rank([ones(size(ind,1),1),log(cli_no_mat(ind,2)),log(cli_no_mat(ind,2)).^2]) == 3 %make sure we have full rank
    [b,~,r] = regress(log(spc(ind)),[ones(size(ind,1),1),log(cli_no_mat(ind,2)),log(cli_no_mat(ind,2)).^2]);              
    rMSE = sqrt(sum(r.^2)/size(ind,1));
    mreg = [b(2);b(3);rMSE];
else
    mreg = [10;10;10];
    display('ERROR: Nothing to regress moms.m sales-per-client regression');
end
if isnan(mreg) ~= [0;0;0]
    mreg = [1;1;2];
    display('ERROR: nans in mreg');
end

%average shipments per client per year
mavship = mean(ship_f_spc(cli_no_mat(:,2)>0)./cli_no_mat(cli_no_mat(:,2)>0,2));

else
    mreg = ones(3,1)*1;
    mavship = ones(1,1)*1;
    mstat = ones(2,1)*1;
    clidist = ones(2,1)*1;
    hazrate = ones(5,1)*1;
    mavex = ones(5,1)*1;
    mnumex = ones(5,1)*1;
    vtran = ones(8,1)*1;
    mexshr = ones(1,1)*1;
    mexreg = ones(3,1)*1;
end
else
    mreg = ones(3,1)*1;
    mavship = ones(1,1)*1;
    mstat = ones(2,1)*1;
    clidist = ones(2,1)*1;
    hazrate = ones(5,1)*1;
    mavex = ones(5,1)*1;
    mnumex = ones(5,1)*1;
    vtran = ones(8,1)*1;
    mexshr = ones(1,1)*1;
    mexreg = ones(3,1)*1;
end
%save test_momsend
end
