% This code snippet reads in initial demand shocks for a single firm, into its 'client hotel' matrix

%create current client count (you need the current client count to know 
%in which column to put the demand shock trajectory)
ind(1:obin,9) = cumsum(ind(1:obin,9)); %home client count
ind(1:obin,10) = cumsum(ind(1:obin,10)); %foreign client count
    
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
