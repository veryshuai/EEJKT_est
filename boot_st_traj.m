function [ind,st,ds,sh,act,breakflag,deathmat] = boot_st_traj(th_ind,mu_f,mu_h,sp_p,lambda_f,lambda_h,c_val_h,c_val_f,mm,ri,rv,rs)

burn    = mm.burn;    %number of burn in periods
esT     = 18 + burn;  % number of ergodic state periods to be simulated
delta   = mm.delta;    %exogenous match death hazard
d       = mm.d;        %exogenous firm death hazard
S       = mm.S;   %number of firms

n_size      = mm.n_size;        %number of matches which are learned from
net_size    = mm.net_size;      %max number of network effects.

Z           = mm.Z;        %grid for other firm productivity
Phi         = mm.Phi;      %grid for own productivity
X_f         = mm.X_f;      %grid for foreign macro shock
X_h         = mm.X_h;      %grid for home macro shock
actual_h    = mm.actual_h; %actual home macro shock indexes (and yrs)
actual_f    = mm.actual_f; %actual foreign macro shock indexes (and yrs)

L_b         = mm.L_b;
L_z         = mm.L_z;
L_f         = mm.L_f;
L_h         = mm.L_h;
erg_pz      = mm.erg_pz;    %stationary distribution of buyer productivities
erg_pp      = mm.erg_pp;    %stationary distribution of seller productivites

%This function does the simulation needed to calculate moments

maxc            = mm.maxc; %maximum number of current clients (follows old program)
max_client_prod = mm.max_client_prod; %maximum changes in demand shock over relationship
mult_match_max  = mm.mult_match_max; %maximum number of matches per exogenous state change interval
mms             = mm.mms; %maximum number of matrix rows (memory issues)

breakflag = 0; %this flag goes to one if there is a maximum matrix violation, and allows us to stop the loop  

match_violation = 0; %this counts number of match per period violations
max_mat_violation = 0; %counts number of matrix size violations
violation = 0; %will count number of demand shock change violations
match_number_violation = 0; %counts number of maximum match number violations
no_more_rands = 0; %counts number of times end of random vector reached 

%% SETTING UP

st = cell(S,1);
ind = cell(S,1);
ds = cell(S,1);
sh = cell(S,1);
act = cell(S,1);
deathmat = cell(S,1);

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

TT = ones(S,1)*esT;

%% EXOGENOUS TRAJECTORIES (MACRO STATES AND SELF PRODUCTIVITY)

obin = ones(S,1)*(esT-burn+2); %current observation index
breakflag = 0;
parfor j = 1:S

    rng(j);
    
    match_violation = 0; %parfor doesn't like keeping track of this
    max_mat_violation = 0; %parfor doesn't like keeping track of this either.
    display(j);
    %if breakflag == 0
    %preallocate
    
    st{j} = repmat(zeros(1,12),mms,1);
    ind{j} = repmat([0,0,0,0,-1,-1,-1,-1,zeros(1,2),-1,-1],mms,1);
    deathmat{j} = zeros(mms,1);
    
    %put in "actual" X_f's (foreign macro shocks)
    for k = 1:size(actual_f,1)-2
        ind{j}(1+k,1) = burn+k;
        ind{j}(1+k,4)  = actual_f(k+2,2);
    end

    %put in "actual" X_h's (home macro shocks)
    for k = 1:size(actual_h,1)-2
        ind{j}(1+k,1) = burn+k;
        ind{j}(1+k,3)  = actual_h(k+2,2);
    end 

    %get Phi's (self productivities)    
    ind{j}(obin(j),2) = find(rand<cum_erg_pp,1,'first');%initial productivity drawn out of ergodic distribution 
    lag = 1;
    time = 0;
    while time < TT
        time = time + exprnd(1/d);
        [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
        [~,death] = find([0,cum_sp_p(j,:)<time],1,'last'); 
        if death ~= 1
            for k = lag:death-1
                ind{j}(obin(j),:) = [cum_sp_p(j,k),...
                ind{j}(obin(j)-1,2)+1-2*(rand<.5*(1+(ind{j}(obin(j)-1,2)-Phi_size-1)/Phi_size)),...
                0,0,-1,-1,-1,-1,zeros(1,2),-1,-1];
                [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end        
            end
            lag = death;
        end
        ind{j}(obin(j),:) = [time,find(rand<cum_erg_pp,1,'first'),0,0,0,0,0,0,zeros(1,2),-1,-1];
        deathmat{j}(obin(j),1) = 1;
    end

    %put in burn period macro shocks (backwards!)

    %foreign
    next = actual_f(3,2);
    gap = exprnd(1/L_f);
    while gap<burn+1
        ind{j}(obin(j),1) = burn+1-gap;
        ind{j}(obin(j),4) = next+1-2*(rand<.5*(1+(next-x_size-1)/x_size));
        next = ind{j}(obin(j),4);
        [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
        gap = gap + exprnd(1/L_f);
    end
    ind{j}(obin(j),1) = 2e-12; %put in correct starting value
    ind{j}(obin(j),4) = next; 
    [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; end

    %home   
    next = actual_h(3,2);
    g = 1;
    gap = exprnd(1/L_h);
    while gap<burn+1
        ind{j}(obin(j),1) = burn+1-gap;
        ind{j}(obin(j),3) = next+1-2*(rand<.5*(1+(next-x_size-1)/x_size));
        next = ind{j}(obin(j),3);
        [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
        gap = gap + exprnd(1/L_h);
    end
    ind{j}(obin(j),1) = 1e-12; %put in correct starting value
    ind{j}(obin(j),3) = next; 
    [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; end

    %% MATCHES


    %resort and fill in
    %first part eliminates times greater than TT, and resets obin(j)
    [obin(j),~] = find(ind{j}(1:obin(j),1)>0,1,'last');
    deathmat{j}(obin(j)+1:end) = 0;
    ind{j}(obin(j)+1:end,:) = repmat([0,0,0,0,-1,-1,-1,-1,...        
    zeros(1,2),-1,-1],size(ind{j}(obin(j)+1:end,:),1),1);
    [ind{j}(1:obin(j),:),I] = sortrows(ind{j}(1:obin(j),:),1);
    deathmat{j}(1:obin(j),1) = deathmat{j}(I,1); 
    for t = 2:4
        [loop_ind,~] = find(ind{j}(2:obin(j),t)==0);
        for k = 1:size(loop_ind) 
            temp = ind{j}(loop_ind(k),t);
            ind{j}(loop_ind(k)+1,t) = temp;
        end
    end

    obin_fix = zeros(S,1);
    %get foreign matches
    obin_fix(j) = obin(j);
    deathind = find(deathmat{j} == 1); 
    lag = 5;
    deathind_ext = [deathind;obin_fix(j)]; 
    for n = 1:size(deathind_ext,1)
           %initialize
           m=0;%match
           s=0;%success
           m_obs = 0; %observed match (able to learn only from observed matches)
           s_obs = 0; %observed success  
           for k = lag:deathind_ext(n)
                %display([j,k,n,ind{j}(k-1,2),ind{j}(k-1,4)]);
                spell = exprnd(1/lambda_f{s_obs+1,m_obs+1,th_ind(j,1),min(s,net_size)+1}...
                (ind{j}(k-1,2),ind{j}(k-1,4))); %time before match
                gap = ind{j}(k,1)-ind{j}(k-1,1); %time between state changes
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
                                rel_time = min(exprnd(1/delta),ind{j}(deathind_ext(n),1)-ind{j}(k-1,1)-cum_spell-1e-12); %exogenous match separation time
                                [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                                ind{j}(obin(j),:) = [ind{j}(k-1,1)+cum_spell,0,0,0,-1,-1,m,s,0,1,0,rel_time];
                                [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                                ind{j}(obin(j),:) = [ind{j}(k-1,1)+cum_spell+rel_time,0,0,0,-1,-1,-1,-1,0,-1,0,0];
                        else
                            [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                            ind{j}(obin(j),:) = [ind{j}(k-1,1)+cum_spell,0,0,0,-1,-1,m,s,0,0,0,0];    
                        end
                        gap = gap - spell; %shrink the gap appropriately 
                        spell = exprnd(1/lambda_f{s_obs+1,m_obs+1,th_ind(j,1),min(s,net_size)+1}(ind{j}(k-1,2),ind{j}(k-1,4))); %swing again
                        cum_spell = spell + cum_spell;
                    end
                end
                lag = deathind_ext(n)+1;    
            end
    end

    %get home matches
    deathind = find(deathmat{j} == 1); 
    lag = 5; %how many rows down the matrix to start (there are some garbage rows at the beginning)
    deathind_ext = [deathind;obin_fix(j)]; 
    for n = 1:size(deathind_ext,1)
        m=0;%match
        s=0;%success
        for k = lag:deathind_ext(n)
            %display([j,k,n,lag,deathind_ext(n),ind{j}(k-1,2),ind{j}(k-1,4),th_ind(j,1),th_ind(j,2),min(s,net_size),ind{j}(k-1,2),ind{j}(k-1,3)]);
            spell = exprnd(1/lambda_h{th_ind(j,1),th_ind(j,2),min(s,net_size)+1}(ind{j}(k-1,2),ind{j}(k-1,3)));
            gap = ind{j}(k,1)-ind{j}(k-1,1);
            if spell < gap             
                p = 0;
                cum_spell = spell;
                while spell < gap && p < mult_match_max
                    p = p+1;
                    [~,match_violation] = step(p,mult_match_max,match_violation);
                    m = m+1;
                    if rand<mu_h(j)
                        s = s + 1;
                        [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                        rel_time = min(exprnd(1/delta),ind{j}(deathind_ext(n),1)-ind{j}(k-1,1)-cum_spell-1e-12);
                        [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                        ind{j}(obin(j),:) = [ind{j}(k-1,1)+cum_spell,0,0,0,m,s,-1,-1,1,0,rel_time,0];
                        [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                        ind{j}(obin(j),:) = [ind{j}(k-1,1)+cum_spell+rel_time,0,0,0,-1,-1,-1,-1,-1,0,0,0];
                    else
                        [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                        [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                        ind{j}(obin(j),:) = [ind{j}(k-1,1)+cum_spell,0,0,0,m,s,-1,-1,0,0,0,0];
                    end  
                    gap = gap - spell;
                    spell = exprnd(1/lambda_h{th_ind(j,1),th_ind(j,2),min(s,net_size)+1}(ind{j}(k-1,2),ind{j}(k-1,3)));
                    cum_spell = spell + cum_spell;
                end
            end
            lag = deathind_ext(n)+1;
        end
    end


    %% Client Trajectories

    %This section creates demand shock trajectories for each successful match
    %recorded in the last section.

    ds{j} = zeros(mms,maxc*2+2);

    [ind{j}(1:obin(j),:),I] = sortrows(ind{j}(1:obin(j),:),1);
    deathmat{j}(1:obin(j),:) = deathmat{j}(I,:);

    %create current client count (you need the current client count to know 
    %in which column to put the demand shock trajectory)
    ind{j}(1:obin(j),9) = cumsum(ind{j}(1:obin(j),9)); %home client count
    ind{j}(1:obin(j),10) = cumsum(ind{j}(1:obin(j),10)); %foreign client count

    %Demand Shock Paths
    g=0; %index for random vector used to calculate client productivities
    occupied = zeros(size(ind{j},1),maxc+1); %0 if slot is unoccupied

    %foreign
    obin_fix = obin(j);
    for k = 1:obin_fix
        if ind{j}(k,12) > 0 && ind{j}(k,10) <= maxc %check that we still have room for more trajectories
            if ind{j}(k,10) == maxc breakflag = 1; match_number_violation = match_number_violation + 1; end %record match number violation
            [~,slot] = find(occupied(k,:)< 1,1,'first'); %finds the next unoccupied slot
            slot = maxc + slot;
            [occ_ind,~] = find(ind{j}(k+1:end,1)>=ind{j}(k,1)+ind{j}(k,12),1,'first');
            if isempty(occ_ind) == 1;
                occ_ind = obin_fix-k;
            end
            occupied(k:k+occ_ind,slot-maxc) = 1; %block off the slot until exogenous match death
            if isempty(slot) == 0
                ds{j}(k,end) = slot(1,1); %record slot number
            else
                slot = maxc;
                display('WARNING: Funny business in demand shock calculation');
            end
            [ds{j}(k,slot),~] = find(rand<cum_erg_pz,1,'first'); %get initial demand shock from ergodic dist.
            cum_spell = exprnd(1/L_z); %time before next demand shock change 

            p = 0;%this keeps track of number of demand shock changes
            while cum_spell < ind{j}(k,12) && p<=max_client_prod %until we reach exogenous separation time, keep the ds changes coming.
                if p == max_client_prod violation = violation + 1; end %record max demand shock change violations
                p = p+1;
                [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                ind{j}(obin(j),:) = [ind{j}(k,1)+cum_spell,0,0,0,-1,-1,-1,-1,-1,-1,0,0]; %make a "blank" row
                ds{j}(obin(j),:) = zeros(1,maxc*2+2);
                if p == 1
                    new_st = ds{j}(k,slot)+1-2*(rand<.5*(1+(ds{j}(k,slot)-z_size-1)/z_size)); 
                else
                    new_st = ds{j}(obin(j)-1,slot)+1-2*(rand<.5*(1+(ds{j}(obin(j)-1,slot)-z_size-1)/z_size));
                end
                ds{j}(obin(j),slot) =  new_st; %record the new state in the ds matrix
                cum_spell = cum_spell + exprnd(1/L_z);
            end
        end
    end

    %home
    occupied = zeros(size(ind{j},1),maxc+1); %reset to zero
    g=0; %index for random vector used to calculate client productivities    
    for k = 1:obin_fix
        if ind{j}(k,11) > 0 && ind{j}(k,9) <= maxc
            if ind{j}(k,9) == maxc breakflag = 1; match_number_violation = match_number_violation + 1; end
            [~,slot] = find(occupied(k,:)< 1,1,'first');
            [occ_ind,~] = find(ind{j}(k+1:end,1)>=ind{j}(k,1)+ind{j}(k,11),1,'first');
            if isempty(occ_ind) == 1;
                occ_ind = obin_fix-k;
            end
            occupied(k:k+occ_ind,slot) = 1;
            if isempty(slot) == 0
                ds{j}(k,end-1) = slot(1,1);
            else
                slot = maxc;
                display('WARNING: Funny business in demand shock calculation');
            end
            ds{j}(k,end-1) = slot(1,1);
            ds{j}(k,slot) = find(rand<cum_erg_pz,1,'first'); %get initial prod from ergodic dist.
            cum_spell = exprnd(1/L_z);  
            p = 0;%this keeps track of number of demand shock changes
            while cum_spell < ind{j}(k,11) && p<=max_client_prod
                if p == max_client_prod violation = violation + 1; end
                p = p+1;
                [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                ind{j}(obin(j),:) = [ind{j}(k,1)+cum_spell,0,0,0,-1,-1,-1,-1,-1,-1,0,0];
                ds{j}(obin(j),:) = zeros(1,maxc*2+2);
                if p == 1
                    new_st = ds{j}(k,slot)+1-2*(rand<.5*(1+(ds{j}(k,slot)-z_size-1)/z_size));
                else
                    new_st = ds{j}(obin(j)-1,slot)+1-2*(rand<.5*(1+(ds{j}(obin(j)-1,slot)-z_size-1)/z_size));
                end
                ds{j}(obin(j),slot) =  new_st;
                cum_spell = cum_spell + exprnd(1/L_z);
            end
        end
    end


    %now resort and fill in
    [ind{j}(1:obin(j),:),I] = sortrows(ind{j}(1:obin(j),:),1);
    ds{j}(1:obin(j),:) = ds{j}(I,:);
    deathmat{j}(1:obin(j),:) = deathmat{j}(I,:);
    ind{j}(1,5:12)=zeros(1,8);
    for t = 2:4
        [loop_ind,~] = find(ind{j}(2:obin(j),t)==0);
        for k = 1:size(loop_ind,1)
            temp = ind{j}(loop_ind(k),t);
            ind{j}(loop_ind(k)+1,t) = temp;
        end
    end
    for t = 5:10
        [loop_ind,~] = find(ind{j}(2:obin(j),t)==-1);
        for k = 1:size(loop_ind,1)
            temp = ind{j}(loop_ind(k),t);
            ind{j}(loop_ind(k)+1,t) = temp;        
        end
    end
    for t = 11:12
        [loop_ind,~] = find(ind{j}(2:obin(j),t)==-1);
        for k = 1:size(loop_ind,1)
            ind{j}(loop_ind(k)+1,t) = 0;  
        end
    end
    %% fill in demand shocks
    %so far we have just recorded the demand shock in the period in which they change.  To
    %calculate profits, we need to fill in the ds matrix.
    %home
    for t = 1:maxc
         [loop_ind,~] = find(ds{j}(:,end-1) == t);
         if isempty(loop_ind) == 0 
             for k = loop_ind'
                rel_time = ind{j}(k,11); %time before exogenous separation
                p = k+1;
                while p<size(ind{j},1) && (ind{j}(p,1)-ind{j}(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
                    if ds{j}(p,t) == 0; 
                             temp = ds{j}(p-1,t);
                             ds{j}(p,t) = temp;    
                    end
                    p = p+1;
                end
            end
        end
    end
    %foreign
    for t = 1:maxc
            [loop_ind,~] = find(ds{j}(:,end) == maxc+t);
            if isempty(loop_ind) == 0
                for k = loop_ind'
                    rel_time = ind{j}(k,12); %time before exogenous separation
                    p = k+1;
                    while p<size(ind{j},1) && (ind{j}(p,1)-ind{j}(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
                        if ds{j}(p,maxc+t) == 0; 
                            temp = ds{j}(p-1,maxc+t); 
                            ds{j}(p,maxc+t)=temp;    
                    end
                    p = p+1;
                end
            end
        end
    end

    %clear drw9 drw10;

    %% Add Shipments
    sh{j} = zeros(mms,2*maxc+2);

    %home
    g = 1;
    for t = 1:maxc
             [loop_ind,~] = find(ds{j}(:,end-1) == t);
             if isempty(loop_ind) == 0 
                 for k = loop_ind'
                     rel_time = ind{j}(k,11); %time before exogenous separation
                     sh{j}(k,t) = ind{j}(k,5); %identity of shipment (to which matched number firm)
                     tot_time = 0; %time spent up to last shipment 
                     while tot_time < rel_time
                         next_ship_time = exprnd(1/L_b);
                         if tot_time+next_ship_time<rel_time
                             [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                             ind{j}(obin(j),:) = [ind{j}(k,1)+tot_time+next_ship_time,0,0,0,-1,-1,-1,-1,-1,-1,0,0];
                             ds{j}(obin(j),:) = [-1*ones(1,2*maxc),0,0];
                             sh{j}(obin(j),:) = zeros(1,2*maxc+2);
                             sh{j}(obin(j),t) = sh{j}(k,t);
                        end
                        tot_time = tot_time + next_ship_time;
                    end
                end
            end
    end
    %foreign
    g = 1;
    for t = 1:maxc
             [loop_ind,~] = find(ds{j}(:,end) == maxc+t);
             if isempty(loop_ind) == 0 
                 for k = loop_ind'
                     rel_time = ind{j}(k,12); %time before exogenous separation
                     sh{j}(k,maxc+t) = ind{j}(k,7);
                     tot_time = 0; %time spent up to last shipment 
                     while tot_time < rel_time
                         next_ship_time = exprnd(1/L_b);
                         if tot_time+next_ship_time<rel_time
                             [obin(j),max_mat_violation,vio] = step(obin(j),mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                             ind{j}(obin(j),:) = [ind{j}(k,1)+tot_time+next_ship_time,0,0,0,-1,-1,-1,-1,-1,-1,0,0];
                             ds{j}(obin(j),:) = [-1*ones(1,2*maxc),0,0];
                             sh{j}(obin(j),:) = zeros(1,2*maxc+2);
                             sh{j}(obin(j),maxc+t) = sh{j}(k,maxc+t);
                        end
                        tot_time = tot_time + next_ship_time;
                    end
                end
            end
    end

    %sort and fill in
    [ind{j}(1:obin(j),:),I] = sortrows(ind{j}(1:obin(j),:),1);
    ds{j}(1:obin(j),:) = ds{j}(I,:);
    sh{j}(1:obin(j),:) = sh{j}(I,:);
    deathmat{j}(1:obin(j),:) = deathmat{j}(I,:);

    for t = 2:4
            [loop_ind,~] = find(ind{j}(2:obin(j),t)==0);
            for k = 1:size(loop_ind,1)
                temp = ind{j}(loop_ind(k),t);
                ind{j}(loop_ind(k)+1,t) = temp;
            end
        end
        for t = 5:10
            [loop_ind,~] = find(ind{j}(2:obin(j),t)==-1);
            for k = 1:size(loop_ind,1)
                temp = ind{j}(loop_ind(k),t);
                ind{j}(loop_ind(k)+1,t) = temp;        
            end
        end
        for t = 1:2*maxc
            [loop_ind,~] = find(ds{j}(2:obin(j),t)==-1);
            for k = 1:size(loop_ind,1)
                ds{j}(loop_ind(k)+1,t) = ds{j}(loop_ind(k),t);
            end
        end

        %% Endogenous separation

        %in this section, we check the continuation value for a particular state to
        %see if it is zero.  If the continuation value is zero, the relationship is
        %immediately endogenously dropped.

        %home    
        deathind = find(deathmat{j} == 1);
        deathind_ex = [deathind;obin(j)];
        lag = 1;
        for n = 1:size(deathind_ex)
            for t = 1:maxc
                [loop_ind,~] = find(ds{j}(lag:deathind_ex(n),end-1) == t);
                loop_ind = loop_ind + lag - 1;
                if isempty(loop_ind) == 0 
                    for k = loop_ind'
                        rel_time = ind{j}(k,11); %time before exogenous separation
                        p = k+1;
                        dropped = 0; %dummy changes to one if relationship is dropped
                        while p<deathind_ex(n) && (ind{j}(p,1)-ind{j}(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
                            if dropped == 0 && sh{j}(p-1,t) ~= 0; %sale is made
                                if c_val_h{ds{j}(p-1,t)}(ind{j}(p-1,2),ind{j}(p-1,3))==0 %check for endogenous separation in the last period
                                    sh{j}(p,t) = 0; %kill any sales in current period  
                                    ind{j}(p,9) = ind{j}(p,9)-1; %reduce current clients by one
                                    dropped = 1;
                            end
                        elseif dropped == 1   
                            sh{j}(p,t) = 0; %kill further sales if relationship is over
                            ind{j}(p,9) = ind{j}(p,9)-1; %reduce current clients by one
                        end
                        p = p+1;
                    end
                    ind{j}(k,11) = ind{j}(find(sh{j}(1:deathind_ex(n),t)==ind{j}(k,5),1,'last'),1)-ind{j}(k,1)+1e-12;
                end
            end
        end
        sh{j}(lag:deathind_ex(n),2*maxc+1) = sum(sh{j}(lag:deathind_ex(n),1:maxc),2); %fill in penultimate column (shipment identities)
        sh{j}(lag:deathind_ex(n),1:maxc) = sh{j}(lag:deathind_ex(n),1:maxc)>0; %make shipments binary (1 for a shipment)

        %foreign
        for t = 1:maxc
                 [loop_ind,~] = find(ds{j}(lag:deathind_ex(n),end) == maxc+t);
                 loop_ind = loop_ind + lag - 1;
                 if isempty(loop_ind) == 0 
                     for k = loop_ind'
                         rel_time = ind{j}(k,12); %time before exogenous separation
                         p = k+1;
                         dropped = 0; %dummy changes to one if relationship is dropped
                         while p<deathind_ex(n) && (ind{j}(p,1)-ind{j}(k,1)<rel_time-1e-12) %note: I had to take a little bit off here to deal with rounding error
                             if dropped == 0 && sh{j}(p-1,maxc+t) ~= 0; %sale is made
                                 if c_val_f{ds{j}(p-1,maxc+t)}(ind{j}(p-1,2),ind{j}(p-1,4))==0 %check for endogenous separation in the last period
                                     sh{j}(p,maxc+t) = 0; %kill any sales in current period  
                                     ind{j}(p,10) = ind{j}(p,10)-1; %reduce current clients by one
                                     dropped = 1;
                            end
                        elseif dropped == 1   
                            sh{j}(p,maxc+t) = 0; %kill further sales if relationship is over
                            ind{j}(p,10) = ind{j}(p,10)-1; %reduce current clients by one
                        end
                        p = p+1;
                    end
                    %display([j,k,t]);
                    %display(ind{j}(find(sh{j}(:,t)==ind{j}(k,5),1,'last'),1)-ind{j}(k,1)+1e-12);
                    ind{j}(k,12) =  ind{j}(find(sh{j}(1:deathind_ex(n),maxc+t)==ind{j}(k,7),1,'last'),1)-ind{j}(k,1)+1e-12;
                end
            end
        end
        sh{j}(lag:deathind_ex(n),2*maxc+2) = sum(sh{j}(lag:deathind_ex(n),maxc+1:2*maxc),2); %fill in last column (shipment identities)
        sh{j}(lag:deathind_ex(n),maxc+1:2*maxc) = sh{j}(lag:deathind_ex(n),maxc+1:2*maxc)>0; %make shipments binary (1 for a shipment)

        lag = deathind_ex(n)+1;
    end

    %% Active years (needed to calculate some moments)
    %active years for each client relationship (NOTE: if a relationship is ongoing,
    %we count the client EVEN IF NO SALES OCCUR during the particular
    %year.)
    act{j} = zeros(mms,2);
    [loop_ind,~] = find(ind{j}(1:obin(j),11)~=0);
    for k = loop_ind
            act{j}(k,1) = max(ceil(ind{j}(k,1)+ind{j}(k,11)),burn)-max(floor(ind{j}(k,1)),burn); %number of active years home
        end
        [loop_ind,~] = find(ind{j}(1:obin(j),12)~=0);
        for k = loop_ind
            act{j}(k,2) = max(ceil(ind{j}(k,1)+ind{j}(k,12)),burn)-max(floor(ind{j}(k,1)),burn); %number of active years foreign
        end

        %% Eliminate unused observations
        [loop_ind,~] = find(ind{j}(:,1) <= TT(j) & ind{j}(:,1)~= 0);
        st{j} = sparse(st{j}(loop_ind,:));
        ds{j} = sparse(ds{j}(loop_ind,:));
        sh{j} = sparse(sh{j}(loop_ind,:));
        act{j} = sparse(act{j}(loop_ind,:));
        ind{j} = sparse(ind{j}(loop_ind,:));
        deathmat{j} = sparse(deathmat{j}(loop_ind,:));

        %% create values array 
        Z_big = [-inf;Z];%tricky way to get the right value for Z          
        if size(ind{j},1)>3 %one rowed ind{j} crashes the program for some reason (situation if a firm dies before burn period ends)
            st{j} = zeros(size(ind{j}));
            st{j}(4:end,:) = [ind{j}(4:end,1),exp(Phi(ind{j}(4:end,2))),exp(X_h(ind{j}(4:end,3))),exp(X_f(ind{j}(4:end,4))),ind{j}(4:end,5:12)];
     else
         st{j} = ind{j};
     end
     for k = 1:size(ds{j},2)-2 %for loop inserted to avoid memory errors
         temp = exp(Z_big(ds{j}(:,k)+1));
         ds{j}(:,k) = temp;
     end
end

% display(['A total of ', num2str(violation),' productivity change violations']);
% display(['A total of ', num2str(match_violation),' match per period violations']);
% display(['A total of ', num2str(match_number_violation),' maximum match violations']);
% display(['A total of ', num2str(max_mat_violation),' matrix size violations']);
% display(['Ran out of rands ', num2str(no_more_rands),' times']);

display(['At end of st_traj, ri is', num2str(ri), 'and the final random draw is ', num2str(rand),'.'])

if breakflag == 1
    display('WARNING: Broke out of loop! Results not reliable.')
end
end
