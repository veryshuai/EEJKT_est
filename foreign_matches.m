% This script snippet creates foreign matches for a single firm

% fix the 'first' index
obin_fix = zeros(S,1);     

%%%%%%%%%%%%%%%%%%%%%GET FOREIGN MATCHES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        if ind(k-1,1) < shock_year % before macro policy change
            lam = lambda_f_orig(s_obs+1,m_obs+1,1,min(s,net_size)+1,ind(k-1,2),ind(k-1,4));
            exp_inv_temp = log(rand)/-lam;
        else
            lam = lambda_f_new(s_obs+1,m_obs+1,1,min(s,net_size)+1,ind(k-1,2),ind(k-1,4));
            exp_inv_temp = log(rand)/-lam;
        end
        cost_vec(k-1,1) = cost(lam,min(s,net_size)+1);
        spell = exp_inv_temp; %time before match
        gap = ind(k,1)-ind(k-1,1); %time between state changes
        if spell < gap
            p=0; %mult match counter
            cum_spell = spell; %cumulative spell since last exogenous state change
            while spell < gap && p < mult_match_max
               p = p+1;
                [~,match_violation] = step(p,mult_match_max,match_violation);
                m = m+1;                
                if rand<succ_prob(obin_fix(j),1) %check for success
                    s = s + 1;
                    if m<=n_size m_obs = m; s_obs = s; end
                    exp_inv_temp = log(rand)/-delta;
                    rel_time = min(exp_inv_temp,ind(deathind(n),1)-ind(k-1,1)-cum_spell-1e-12); %exogenous match separation time
                    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                    ind(obin,:) = [ind(k-1,1)+cum_spell,0,0,0,-1,-1,m,s,0,1,0,rel_time];
                    cost_vec(obin,2) = F;
                    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                    ind(obin,:) = [ind(k-1,1)+cum_spell+rel_time,0,0,0,-1,-1,-1,-1,0,-1,0,0];
                else
                    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                    ind(obin,:) = [ind(k-1,1)+cum_spell,0,0,0,-1,-1,m,s,0,0,0,0];    
                end
                gap = gap - spell; %shrink the gap appropriately 
                if ind(k-1,1) + cum_spell < shock_year % before macro policy change
                    lam = lambda_f_orig(s_obs+1,m_obs+1,1,min(s,net_size)+1,ind(k-1,2),ind(k-1,4));
                    exp_inv_temp = log(rand)/-lam; %swing again
                else
                    lam=lambda_f_new(s_obs+1,m_obs+1,1,min(s,net_size)+1,ind(k-1,2),ind(k-1,4));
                    exp_inv_temp = log(rand)/-lam; %swing again
                end
                cost_vec(obin,1) = cost(lam,min(s,net_size)+1);
                spell = exp_inv_temp;
                cum_spell = spell + cum_spell;
            end
        end
    lag = deathind(n)+1;
    end
end
