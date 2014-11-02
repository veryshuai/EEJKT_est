function [sh_val_h,sh_val_f,ind,deathmat,ds,sh,act,breakflag,max_mat_violation,match_number_violation,prod_init] = singlefirm(j,max_mat_violation,match_violation,violation,match_number_violation,no_more_rands,breakflag,x_size,Phi_size,z_size,TT,cum_erg_pz,cum_erg_pp,cum_sp_p,th_ind,mu_h,mu_f,sp_p,c_val_h_orig,c_val_f_orig,c_val_h_new,c_val_f_new,lambda_f_orig,lambda_h_orig,lambda_f_new,lambda_h_new,burn,delta,d,S,n_size,net_size,Z,Phi,X_f,X_h,actual_h,actual_f,L_b,L_z,L_f,L_h,erg_pz,erg_pp,maxc,max_client_prod,mult_match_max,mms,scale_h,scale_f,de,agg_shocks)
% This function generates the state matrices for a single firm, to be optionally
% compilied into mex.

% repeatedly reseed in the sequel:
seed = randi(1e7,6,1);

Z_big = [-inf,Z'];

coder.extrinsic('display');
exp_inv_temp = 0;
scalar_temp = 0;

% display(j);
obin = 20; %current observation index
st       = zeros(mms,12);
ind      = zeros(mms,12);
deathmat = zeros(mms,1);
ds       = zeros(mms,maxc*2+2);
sh       = zeros(mms,maxc*2+2);
act      = zeros(mms,2);
sh_val_h = zeros(mms,3);
sh_val_f = zeros(mms,3);

if breakflag == 0

    %preallocate
    ind(:,:) = repmat([0,0,0,0,-1,-1,-1,-1,zeros(1,2),-1,-1],mms,1);
    deathmat(:) = zeros(mms,1);
    
    %put in "actual" X_f's (foreign macro shocks)
    for k = 1:size(actual_f,1)-2
        ind(1+k,1) = burn+k;
        ind(1+k,4)  = actual_f(k+2,2);
    end

    %put in "actual" X_h's (home macro shocks)
    for k = 1:size(actual_h,1)-2
        ind(1+k,1) = burn+k;
        ind(1+k,3)  = actual_h(k+2,2);
    end 

    %get Phi's (self productivities)
    prod_init = find(rand<cum_erg_pp,1,'first');%initial productivity drawn out of ergodic distribution 
    scalar_prod_init = prod_init(1,1);
    ind(obin,2) = scalar_prod_init;
    lag = 1;
    time = 0;
    while time < TT
        exp_inv_temp = log(rand)/-d;
        time = time + exp_inv_temp; %time of death
        [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
        death = find([0,cum_sp_p(j,lag:end)]'<time,1,'last'); %index of death
        scalar_temp = death(1,1)-1; 
        if scalar_temp ~= 0
            for k = lag:lag+scalar_temp-1
                ind(obin,:) = [cum_sp_p(j,k),...
                    ind(obin-1,2)+1-2*(rand<.5*(1+(ind(obin-1,2)-Phi_size-1)/Phi_size)),...
                    0,0,-1,-1,-1,-1,zeros(1,2),-1,-1];
                [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
            end
            lag = scalar_temp+1;
        end
        ind(obin,:) = [time,find(rand<cum_erg_pp,1,'first'),0,0,0,0,0,0,zeros(1,2),-1,-1];
        deathmat(obin,1) = 1; %new firm
    end

%put in burn period macro shocks (backwards!)
% NEED TO FIX THIS SO THAT ALL FIRMS FACE SAME SHOCKS!
%foreign
k = 1;
[obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; end
next = actual_f(3,2);
exp_inv_temp = log(agg_shocks(k))/-L_f;
k = k+1;
gap = exp_inv_temp;
    while gap<burn+1
        ind(obin,1) = burn+1-gap;
        ind(obin,4) = next+1-2*(agg_shocks(k)<.5*(1+(next-x_size-1)/x_size));
        k = k+1;
        next = ind(obin,4);
        [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
        exp_inv_temp = log(agg_shocks(k))/-L_f;
        k = k+1;
        gap = gap + exp_inv_temp;
    end
ind(obin,1) = 2e-12; %put in correct starting value
ind(obin,4) = next; 
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

%% MATCHES

%resort and fill in
    %first part eliminates times greater than TT, and resets obin
    obin = find((ind(1:obin,1))>0,1,'last');
    obin = obin(1,1);
    deathmat(obin+1:end) = 0;
    scalar_temp = numel(obin+1:mms);
    ind(obin+1:end,:) = repmat([0,0,0,0,-1,-1,-1,-1,zeros(1,2),-1,-1],scalar_temp,1);
    [ind(1:obin,:),I] = sortrows((ind(1:obin,:)),1);
    deathmat(1:obin,1) = deathmat(I,1); 
    for t = 2:4
        [loop_ind,~] = find(ind(2:obin,t)==0);
        for k = 1:size(loop_ind,1)
            temp = ind(loop_ind(k),t);
            ind(loop_ind(k)+1,t) = temp;
        end
    end   

rng(seed(1));

obin_fix = zeros(S,1);
%get foreign matches
        obin_fix(j) = obin;
        deathind = find(deathmat(:) == 1); 
        lag = 5;
        deathind = [deathind;obin_fix(j)]; 
        for n = 1:size(deathind,1)
            %initialize
                m=0;%match
                s=0;%success
                m_obs = 0; %observed match (able to learn only from observed matches)
                s_obs = 0; %observed success  
            for k = lag:deathind(n)
                if ind(k-1,1) < TT - 9 % before macro policy change
                    exp_inv_temp = log(rand)/-lambda_f_orig(s_obs+1,m_obs+1,th_ind(j,1),min(s,net_size)+1,ind(k-1,2),ind(k-1,4));
                else
                    exp_inv_temp = log(rand)/-lambda_f_new(s_obs+1,m_obs+1,th_ind(j,1),min(s,net_size)+1,ind(k-1,2),ind(k-1,4));
                end
                spell = exp_inv_temp; %time before match
                gap = ind(k,1)-ind(k-1,1); %time between state changes
                if spell < gap
                    p=0; %mult match counter
                    cum_spell = spell; %cumulative spell since last exogenous state change
                    while spell < gap && p < mult_match_max
                        p = p+1;
                        [~,match_violation] = step(p,mult_match_max,match_violation);
                        m = m+1;                
                        if rand<mu_f(j) %check for success
                            s = s + 1;
                            if m<=n_size m_obs = m; s_obs = s; end
                            exp_inv_temp = log(rand)/-delta;
                            rel_time = min(exp_inv_temp,ind(deathind(n),1)-ind(k-1,1)-cum_spell-1e-12); %exogenous match separation time
                            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            ind(obin,:) = [ind(k-1,1)+cum_spell,0,0,0,-1,-1,m,s,0,1,0,rel_time];
                            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            ind(obin,:) = [ind(k-1,1)+cum_spell+rel_time,0,0,0,-1,-1,-1,-1,0,-1,0,0];
                        else
                            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            ind(obin,:) = [ind(k-1,1)+cum_spell,0,0,0,-1,-1,m,s,0,0,0,0];    
                        end
                        gap = gap - spell; %shrink the gap appropriately 
                        if ind(k-1,1) < TT - 9 % before macro policy change
                            exp_inv_temp = log(rand)/-lambda_f_orig(s_obs+1,m_obs+1,th_ind(j,1),min(s,net_size)+1,ind(k-1,2),ind(k-1,4)); %swing again
                        else
                            exp_inv_temp = log(rand)/-lambda_f_new(s_obs+1,m_obs+1,th_ind(j,1),min(s,net_size)+1,ind(k-1,2),ind(k-1,4)); %swing again
                        end
                        spell = exp_inv_temp;
                        cum_spell = spell + cum_spell;
                    end
                end
            lag = deathind(n)+1;    
            end
        end

%reseed
rng(seed(2));

%get home matches
        deathind = find(deathmat(:) == 1); 
        lag = 5; %how many rows down the matrix to start (there are some garbage rows at the beginning)
        deathind = [deathind;obin_fix(j)]; 
        for n = 1:size(deathind,1)
            m=0;%match
            s=0;%success
            for k = lag:deathind(n)
                if ind(k-1,1) < TT-9
                    exp_inv_temp = log(rand)/-lambda_h_orig(th_ind(j,1),th_ind(j,2),min(s,net_size)+1,ind(k-1,2),ind(k-1,3));
                else
                    exp_inv_temp = log(rand)/-lambda_h_new(th_ind(j,1),th_ind(j,2),min(s,net_size)+1,ind(k-1,2),ind(k-1,3));
                end
                spell = exp_inv_temp;
                gap = ind(k,1)-ind(k-1,1);
                if spell < gap             
                    p = 0;
                    cum_spell = spell;
                    while spell < gap && p < mult_match_max
                        p = p+1;
                        [~,match_violation] = step(p,mult_match_max,match_violation);
                        m = m+1;
                        if rand<mu_h(j)
                        s = s + 1;
                            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            exp_inv_temp = log(rand)/-delta;
                            rel_time = min(exp_inv_temp,ind(deathind(n),1)-ind(k-1,1)-cum_spell-1e-12);
                            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            ind(obin,:) = [ind(k-1,1)+cum_spell,0,0,0,m,s,-1,-1,1,0,rel_time,0];
                            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            ind(obin,:) = [ind(k-1,1)+cum_spell+rel_time,0,0,0,-1,-1,-1,-1,-1,0,0,0];
                        else
                            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            ind(obin,:) = [ind(k-1,1)+cum_spell,0,0,0,m,s,-1,-1,0,0,0,0];
                        end  
                        gap = gap - spell;
                        if ind(k-1,1) < TT-9
                            exp_inv_temp = log(rand)/-lambda_h_orig(th_ind(j,1),th_ind(j,2),min(s,net_size)+1,ind(k-1,2),ind(k-1,3));
                        else
                            exp_inv_temp = log(rand)/-lambda_h_new(th_ind(j,1),th_ind(j,2),min(s,net_size)+1,ind(k-1,2),ind(k-1,3));
                        end
                        spell = exp_inv_temp;
                        cum_spell = spell + cum_spell;
                    end
                end
            lag = deathind(n)+1;
            end
        end

    
%% Client Trajectories

%This section creates demand shock trajectories for each successful match
%recorded in the last section.

    %ds(:,:) = (zeros(mms,maxc*2+2));
    
    [ind(1:obin,:),I] = sortrows((ind(1:obin,:)),1);
    deathmat(1:obin,:) = deathmat(I,:);
    
    %create current client count (you need the current client count to know 
    %in which column to put the demand shock trajectory)
    ind(1:obin,9) = cumsum(ind(1:obin,9)); %home client count
    ind(1:obin,10) = cumsum(ind(1:obin,10)); %foreign client count

    %Demand Shock Paths
    g=0; %index for random vector used to calculate client productivities
    occupied = zeros(size(ind,1),maxc+1); %0 if slot is unoccupied

%reseed
rng(seed(3));
    
%foreign
obin_fix = obin;
obin_fix = obin_fix(1,1);
for k = 1:obin_fix
    if ind(k,12) > 0 && ind(k,10) <= maxc %check that we still have room for more trajectories
        if ind(k,10) == maxc breakflag = 1; match_number_violation = match_number_violation + 1; end %record match number violation
        [~,slot] = find(occupied(k,:) < 1,1,'first'); %finds the next unoccupied slot
        slot = maxc + slot;
        occ_ind = find((ind(k+1:end,1))>=ind(k,1)+ind(k,12),1,'first');
        if isempty(occ_ind) == 1;
            occ_ind = obin_fix-k;
        end
        scalar_temp = occ_ind(1,1);
        occupied(k:k+scalar_temp,slot-maxc) = 1; %block off the slot until exogenous match death
        if isempty(slot) == 0
            ds(k,end) = (slot(1,1)); %record slot number
        else
            slot = maxc;
            display('WARNING: Funny business in demand shock calculation');
        end
        ds(k,slot) = (find(rand<cum_erg_pz,1,'first')); %get initial demand shock from ergodic dist.
        exp_inv_temp = log(rand)/-L_z; %time before next demand shock change 
        cum_spell = exp_inv_temp;

        p = 0;%this keeps track of number of demand shock changes
        while cum_spell < ind(k,12) && p<=max_client_prod %until we reach exogenous separation time, keep the ds changes coming.
            if p == max_client_prod violation = violation + 1; end %record max demand shock change violations
            p = p+1;
            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
            ind(obin,:) = [ind(k,1)+cum_spell,0,0,0,-1,-1,-1,-1,-1,-1,0,0]; %make a "blank" row
            ds(obin,:) = (zeros(1,maxc*2+2));
            if p == 1
                scalar_temp = 1-2*(rand<.5*(1+(ds(k,slot)-z_size-1)/z_size));
                scalar_temp = scalar_temp(1,1);
                new_st = ds(k,slot)+scalar_temp; 
            else
                scalar_temp = 1-2*(rand<.5*(1+(ds(obin-1,slot)-z_size-1)/z_size));
                scalar_temp = scalar_temp(1,1);
                new_st = ds(obin-1,slot)+scalar_temp;
            end
            ds(obin,slot) =  (new_st); %record the new state in the ds matrix
            exp_inv_temp = log(rand)/-L_z;
            cum_spell = cum_spell + exp_inv_temp;
        end
    end
end

%reseed
rng(seed(4));

%home
occupied = zeros(size(ind,1),maxc+1); %reset to zero
%g=0; %index for random vector used to calculate client productivities    
 for k = 1:obin_fix
    if ind(k,11) > 0 && ind(k,9) <= maxc
        if ind(k,9) == maxc breakflag = 1; match_number_violation = match_number_violation + 1; end
        [~,slot] = find(occupied(k,:)< 1,1,'first');
        occ_ind = find((ind(k+1:end,1))>=ind(k,1)+ind(k,11),1,'first');
        if isempty(occ_ind) == 1;
            occ_ind = obin_fix-k;
        end
        scalar_temp = occ_ind(1,1);
        occupied(k:k+scalar_temp,slot) = 1;
        if isempty(slot) == 0
            ds(k,end-1) = (slot(1,1));
        else
            slot = maxc;
            display('WARNING: Funny business in demand shock calculation');
        end
        ds(k,end-1) = (slot(1,1));
        ds(k,slot) = (find(rand<cum_erg_pz,1,'first')); %get initial prod from ergodic dist.
        exp_inv_temp = log(rand)/-L_z;  
        cum_spell = exp_inv_temp;
        p = 0;%this keeps track of number of demand shock changes
        while cum_spell < ind(k,11) && p<=max_client_prod
            if p == max_client_prod violation = violation + 1; end
            p = p+1;
            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
            ind(obin,:) = [ind(k,1)+cum_spell,0,0,0,-1,-1,-1,-1,-1,-1,0,0];
            ds(obin,:) = (zeros(1,maxc*2+2));
            if p == 1
                scalar_temp = 1-2*(rand<.5*(1+(ds(k,slot)-z_size-1)/z_size));
                scalar_temp = scalar_temp(1,1);
                new_st = ds(k,slot)+scalar_temp;
            else
                scalar_temp = 1-2*(rand<.5*(1+(ds(obin-1,slot)-z_size-1)/z_size));
                scalar_temp = scalar_temp(1,1);
                new_st = ds(obin-1,slot)+scalar_temp;
            end
            ds(obin,slot) =  (new_st);
            exp_inv_temp = log(rand)/-L_z;
            cum_spell = cum_spell + exp_inv_temp;
        end
    end
end
    
    %now resort and fill in
        [ind(1:obin,:),I] = sortrows((ind(1:obin,:)),1);
        ds(1:obin,:) = (ds(I,:));
        deathmat(1:obin,:) = deathmat(I,:);
        ind(1,5:12)=zeros(1,8);
        for t = 2:4
            [loop_ind,~] = find((ind(2:obin,t)==0));
            if isempty(loop_ind) == 0
                for k = 1:size(loop_ind,1)
                    ind(loop_ind(k)+1,t) = ind(loop_ind(k),t);
                end
            end
        end
        for t = 5:10
            [loop_ind,~] = find((ind(2:obin,t))==-1);
            if isempty(loop_ind) == 0
                for k = 1:size(loop_ind,1)
                    temp = ind(loop_ind(k),t);
                    ind(loop_ind(k)+1,t) = temp;        
                end
            end
        end
        for t = 11:12
            [loop_ind,~] = find((ind(2:obin,t))==-1);
            if isempty(loop_ind) == 0
                for k = 1:size(loop_ind,1)
                    ind(loop_ind(k)+1,t) = 0;  
                end
            end
        end

%% fill in demand shocks
 %so far we have just recorded the demand shock in the period in which they change.  To
 %calculate profits, we need to fill in the ds matrix.
    %home
        for t = 1:maxc
             loop_ind = find((ds(:,end-1)) == t);
              if isempty(loop_ind) == 0 
                for o = 1:size(loop_ind,1)
                    k = loop_ind(o);
                rel_time = ind(k,11); %time before exogenous separation
                p = k+1;
                while p<size(ind,1) && (ind(p,1)-ind(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
                    if ds(p,t) == 0; 
                        temp = ds(p-1,t);
                        ds(p,t) = temp;    
                    end
                    p = p+1;
                end
                end
              end
         end
        %foreign
        for t = 1:maxc
            loop_ind = find((ds(:,end)) == maxc+t);
            if isempty(loop_ind) == 0
                for o = 1:size(loop_ind,1)
                    k = loop_ind(o);
                rel_time = ind(k,12); %time before exogenous separation
                p = k+1;
                while p<size(ind,1) && (ind(p,1)-ind(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
                    if ds(p,maxc+t) == 0; 
                        temp = ds(p-1,maxc+t); 
                        ds(p,maxc+t)=temp;    
                    end
                    p = p+1;
                end
                end
            end
        end

%% Add Shipments

%reseed
rng(seed(5));

%seed barrage, need one for each 't'
seed_barrage = randi(1e7,maxc,1);

%home
    g = 1;
    for t = 1:maxc
             rng(seed_barrage(t));
             loop_ind = find((ds(:,end-1)) == t);
              if isempty(loop_ind) == 0 
                for o = 1:size(loop_ind,1)
                    k = loop_ind(o);
                    rel_time = ind(k,11); %time before exogenous separation
                    sh(k,t) = ind(k,5); %identity of shipment (to which matched number firm)
                    sh_val_h(k,2) = ind(k,5); %copy identity into sh_val
                    sh_val_h(k,3) = 1; %signifies that it is a new client
                    tot_time = 0; %time spent up to last shipment 
                    while tot_time < rel_time
                        exp_inv_temp = log(rand)/-L_b;
                        next_ship_time = exp_inv_temp;
                        if tot_time+next_ship_time<rel_time
                            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            ind(obin,:) = [ind(k,1)+tot_time+next_ship_time,0,0,0,-1,-1,-1,-1,-1,-1,0,0];
                            ds(obin,:) = [-1*ones(1,2*maxc),0,0];
                            sh(obin,:) = zeros(1,2*maxc+2);
                            sh(obin,t) = sh(k,t);
                            sh_val_h(obin,2) = ind(k,5); %copy identity into sh_val
                        end
                        tot_time = tot_time + next_ship_time;
                    end
                end
              end
    end

    %reseed
    rng(seed(6));
    
    %seed barrage, need one for each 't'
    seed_barrage = randi(1e7,maxc,1);
    
    %foreign
    g = 1;
    for t = 1:maxc
            rng(seed_barrage(t));
            loop_ind = find((ds(:,end)) == maxc+t);
              if isempty(loop_ind) == 0 
                for o = 1:size(loop_ind,1)
                    k = loop_ind(o);
                    rel_time = ind(k,12); %time before exogenous separation
                    sh(k,maxc+t) = ind(k,7);
                    sh_val_f(k,2) = ind(k,7); %copy identity into sh_val
                    sh_val_f(k,3) = 1; %signifies that it is a new client
                    tot_time = 0; %time spent up to last shipment 
                    while tot_time < rel_time
                        exp_inv_temp = log(rand)/-L_b;
                        next_ship_time = exp_inv_temp;
                        if tot_time+next_ship_time<rel_time
                            [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            ind(obin,:) = [ind(k,1)+tot_time+next_ship_time,0,0,0,-1,-1,-1,-1,-1,-1,0,0];
                            ds(obin,:) = [-1*ones(1,2*maxc),0,0];
                            sh(obin,:) = zeros(1,2*maxc+2);
                            sh(obin,maxc+t) = sh(k,maxc+t);
                            sh_val_f(obin,2) = ind(k,7); %copy identity into sh_val
                        end
                        tot_time = tot_time + next_ship_time;
                    end
                end
              end
    end
        
%sort and fill in
[ind(1:obin,:),I] = sortrows((ind(1:obin,:)),1);
ds(1:obin,:) = ds(I,:);
sh(1:obin,:) = sh(I,:);
sh_val_h(1:obin,:) = sh_val_h(I,:);
sh_val_f(1:obin,:) = sh_val_f(I,:);
deathmat(1:obin,:) = deathmat(I,:);
    
        for t = 2:4
            [loop_ind,~] = find((ind(2:obin,t))==0);
            if isempty(loop_ind)==0
                for k = 1:size(loop_ind,1)
                    temp = ind(loop_ind(k),t);
                    ind(loop_ind(k)+1,t) = temp;
                end
            end
        end
        for t = 5:10
            [loop_ind,~] = find((ind(2:obin,t))==-1);
            for k = 1:size(loop_ind,1)
                temp = ind(loop_ind(k),t);
                ind(loop_ind(k)+1,t) = temp;        
            end
        end
        for t = 1:2*maxc
            loop_ind = find((ds(2:obin,t))==-1);
            for k = 1:size(loop_ind,1)
                ds(loop_ind(k)+1,t) = ds(loop_ind(k),t);
            end
        end
        loop_ind = find(sh_val_h(:,2) ~= 0);
        for k = 1:size(loop_ind,1)
            cli_ind = find(ds(loop_ind(k),:),1,'first');
            sh_val_h(loop_ind(k),1) = exp(scale_h)*exp(Phi(ind(loop_ind(k),2)))^(de-1)*exp(X_h(ind(loop_ind(k),3)))*exp(Z_big(ds(loop_ind(k),cli_ind)+1)); 
        end
        loop_ind = find(sh_val_f(:,2) ~= 0);
        for k = 1:size(loop_ind,1)
            cli_ind = find(ds(loop_ind(k),:),1,'first');
            sh_val_f(loop_ind(k),1) = exp(scale_f)*exp(Phi(ind(loop_ind(k),2)))^(de-1)*exp(X_f(ind(loop_ind(k),3)))*exp(Z_big(ds(loop_ind(k),cli_ind)+1)); 
        end


%% Endogenous separation

%in this section, we check the continuation value for a particular state to
%see if it is zero.  If the continuation value is zero, the relationship is
%immediately endogenously dropped.

%home    
        deathind_before_end = find(deathmat(:) == 1);
        deathind = [deathind_before_end;obin-1];
        lag = 1;
        for n = 1:size(deathind,1)
            for t = 1:maxc
                loop_ind = find((ds(lag:deathind(n),end-1)) == t);
                loop_ind = loop_ind + lag - 1;
                if isempty(loop_ind) == 0 
                    for o = 1:size(loop_ind,1)
                        k = loop_ind(o);
                    rel_time = ind(k,11); %time before exogenous separation
                    p = k+1;
                    dropped = 0; %dummy changes to one if relationship is dropped
                    while p<deathind(n) && (ind(p,1)-ind(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
                        if dropped == 0 && sh(p-1,t) ~= 0; %sale is made
                            if ind(k-1,1) < TT - 9 % before macro policy change
                                if c_val_h_orig(ds(p-1,t),ind(p-1,2),ind(p-1,3))==0 %check for endogenous separation in the last period
                                    sh(p,t) = 0; %kill any shipments in current period  
                                    ind(p,9) = ind(p,9)-1; %reduce current clients by one
                                    dropped = 1;
                                end
                            else
                                if c_val_h_new(ds(p-1,t),ind(p-1,2),ind(p-1,3))==0 %check for endogenous separation in the last period
                                    sh(p,t) = 0; %kill any shipments in current period  
                                    ind(p,9) = ind(p,9)-1; %reduce current clients by one
                                    dropped = 1;
                                end
                            end
                        elseif dropped == 1   
                            sh(p,t) = 0; %kill further shipments if relationship is over
                            ind(p,9) = ind(p,9)-1; %reduce current clients by one
                        end
                        p = p+1;
                    end
                    temp = find((sh(1:deathind(n),t))==(ind(k,5)),1,'last');
                    temp = temp(1,1);
                    ind(k,11) = ind(temp,1)-ind(k,1)+1e-12;
                    end
                end
            end
            sh(lag:deathind(n),2*maxc+1) = sum(sh(lag:deathind(n),1:maxc),2); %fill in penultimate column (shipment identities)
            sh(lag:deathind(n),1:maxc) = (sh(lag:deathind(n),1:maxc))>0; %make shipments binary (1 for a shipment)
        
            %foreign
             for t = 1:maxc
                 loop_ind = find((ds(lag:deathind(n),end)) == maxc+t);
                 loop_ind = loop_ind + lag - 1;
                  if isempty(loop_ind) == 0 
                    for o = 1:size(loop_ind,1)
                        k = loop_ind(o);
                    rel_time = ind(k,12); %time before exogenous separation
                    p = k+1;
                    dropped = 0; %dummy changes to one if relationship is dropped
                    while p<deathind(n) && (ind(p,1)-ind(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
                        if dropped == 0 && sh(p-1,maxc+t) ~= 0; %sale is made
                            if ind(k-1,1) < TT - 9 % before macro policy change
                                if c_val_f_orig(ds(p-1,maxc+t),ind(p-1,2),ind(p-1,4))==0 %check for endogenous separation in the last period 
                                    sh(p,maxc+t) = 0; %kill any sales in current period  
                                    ind(p,10) = ind(p,10)-1; %reduce current clients by one
                                    dropped = 1;
                                end
                            else
                                if c_val_f_new(ds(p-1,maxc+t),ind(p-1,2),ind(p-1,4))==0 %check for endogenous separation in the last period 
                                    sh(p,maxc+t) = 0; %kill any sales in current period  
                                    ind(p,10) = ind(p,10)-1; %reduce current clients by one
                                    dropped = 1;
                                end
                            end
                        elseif dropped == 1   
                            sh(p,maxc+t) = 0; %kill further sales if relationship is over
                            ind(p,10) = ind(p,10)-1; %reduce current clients by one
                        end
                        p = p+1;
                    end
                    %display([k,t]);
                    %display(ind(find(sh(:,t)==ind(k,5),1,'last'),1)-ind(k,1)+1e-12);
                    temp = find((sh(1:deathind(n),maxc+t))==(ind(k,7)),1,'last');
                    temp = temp(1,1);
                    ind(k,12) =  ind(temp,1)-ind(k,1)+1e-12;
                    end
                  end
              end
            sh(lag:deathind(n),2*maxc+2) = sum(sh(lag:deathind(n),maxc+1:2*maxc),2); %fill in last column (shipment identities)
            sh(lag:deathind(n),maxc+1:2*maxc) = (sh(lag:deathind(n),maxc+1:2*maxc))>0; %make shipments binary (1 for a shipment)
            
        lag = deathind(n)+1;
        end
        
        %now kill the sales
        sh_val_h(:,1:2) = [sh_val_h(:,1).*sum(sh(:,1:maxc),2),sh_val_h(:,2).*sum(sh(:,1:maxc),2)];
        sh_val_f(:,1:2) = [sh_val_f(:,1).*sum(sh(:,maxc+1:2*maxc),2),sh_val_f(:,2).*sum(sh(:,maxc+1:2*maxc),2)];

%% Active years (needed to calculate some moments)
    %active years for each client relationship (NOTE: if a relationship is ongoing,
    %we count the client EVEN IF NO SALES OCCUR during the particular
    %year.)
        loop_ind = find((ind(1:obin,11))~=0);
        if isempty(loop_ind)==0
            for o = 1:size(loop_ind,1)
                k = loop_ind(o);
                act(k,1) = max(ceil(ind(k,1)+ind(k,11)),burn)-max(floor(ind(k,1)),burn); %number of active years home
            end
        end
        loop_ind = find((ind(1:obin,12))~=0);
        if isempty(loop_ind)==0
            for o = 1:size(loop_ind,1)
                k = loop_ind(o);
                act(k,2) = max(ceil(ind(k,1)+ind(k,12)),burn)-max(floor(ind(k,1)),burn); %number of active years foreign
            end
        end
    end
end
