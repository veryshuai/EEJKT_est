% Script snippet reads macro shocks in burn in period into even matrix

%put in "actual" X_f's (foreign macro shocks) at end of period
for k = 1:size(actual_f,1)-2 
    ind(k+1,1) = TT - (size(actual_f,1) - 2) + k; % time of the shock
    ind(k+1,4)  = actual_f(k+2,2); % value of the shock
end

%put in "actual" X_h's (home macro shocks) at end of period
for k = 1:size(actual_h,1)-2
    ind(k+1,1) = TT - (size(actual_h,1) - 2) + k; %time of the shock 
    ind(k+1,3)  = actual_h(k+2,2); % value of the shock
end 

%%%%%%%%%%%%%%%%%%%%%Burn in macro shocks%%%%%%%%%%%%%%%%%%%%%%

%foreign

k = 1; % placeholder for shock values in preallocated vector (so all firms get the same macro shock)
[obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; end 
next = actual_f(3,2); %needs to end in the first actual shock
exp_inv_temp = log(agg_shocks(k))/-L_f; %time since last change
k = k+1; 
gap = exp_inv_temp;
while gap<burn+1 %check to make sure there is still time for another shock
    ind(obin,1) = burn+1-gap; %time of shock
    ind(obin,4) = next+1-2*(agg_shocks(k)<.5*(1+(next-x_size-1)/x_size)); %value of shock
    k = k+1;
    next = ind(obin,4);
    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
    exp_inv_temp = log(agg_shocks(k))/-L_f; %time until next shock
    k = k+1;
    gap = gap + exp_inv_temp;
end
ind(obin,1) = 2e-12; %put in correct starting time
ind(obin,4) = next; %put in correct starting value 
[obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; end

%home
next = actual_h(3,2);
g = 1;
exp_inv_temp = log(agg_shocks(k))/-L_h;
k = k+1;
gap = exp_inv_temp;
while gap<burn+1
    ind(obin,1) = burn+1-gap;
    ind(obin,3) = next+1-2*(agg_shocks(k)<.5*(1+(next-x_size-1)/x_size));
    k = k+1;
    next = ind(obin,3);
    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
    exp_inv_temp = log(agg_shocks(k))/-L_h;
    k = k+1;
    gap = gap + exp_inv_temp;
end
ind(obin,1) = 1e-12; %put in correct starting value
ind(obin,3) = next; 
[obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; end
