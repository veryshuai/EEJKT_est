% This script snippet creates home matches for a single firm

%%%%%%%%%%%%%%%%%%GET HOME MATCHES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deathind = find(deathmat(:) == 1); 
lag = 5; %how many rows down the matrix to start (there are some garbage rows at the beginning)
deathind = [deathind;obin_fix(j)]; 
for n = 1:size(deathind,1)
    m=0;%match
    s=0;%success
    for k = lag:deathind(n)
        if ind(k-1,1) < shock_year
            lam = lambda_h_orig(1,succ_prob(obin_fix(j),4),min(s,net_size)+1,ind(k-1,2),ind(k-1,3));
            exp_inv_temp = log(rand)/-lam;
        else
            lam = lambda_h_new(1,succ_prob(obin_fix(j),4),min(s,net_size)+1,ind(k-1,2),ind(k-1,3));
            exp_inv_temp = log(rand)/-lam;
        end
        cost_vec(k,3) = cost(lam,min(s,net_size)+1);
        spell = exp_inv_temp;
        gap = ind(k,1)-ind(k-1,1);
        if spell < gap             
            p = 0;
            cum_spell = spell;
            while spell < gap && p < mult_match_max
                p = p+1;
                [~,match_violation] = step(p,mult_match_max,match_violation);
                m = m+1;
                if rand<succ_prob(obin_fix(j),2) %check for success
                s = s + 1;
                    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                    exp_inv_temp = log(rand)/-delta;
                    rel_time = min(exp_inv_temp,ind(deathind(n),1)-ind(k-1,1)-cum_spell-1e-12);
                    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                    ind(obin,:) = [ind(k-1,1)+cum_spell,0,0,0,m,s,-1,-1,1,0,rel_time,0];
                    cost_vec(obin,4) = F;
                    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                    ind(obin,:) = [ind(k-1,1)+cum_spell+rel_time,0,0,0,-1,-1,-1,-1,-1,0,0,0];
                else
                    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
                    ind(obin,:) = [ind(k-1,1)+cum_spell,0,0,0,m,s,-1,-1,0,0,0,0];
                end  
                gap = gap - spell;
                if ind(k-1,1) + cum_spell < shock_year
                    lam = lambda_h_orig(1,succ_prob(obin_fix(j),4),min(s,net_size)+1,ind(k-1,2),ind(k-1,3));
                    exp_inv_temp = log(rand)/-lam;
                else
                    lam = lambda_h_new(1,succ_prob(obin_fix(j),4),min(s,net_size)+1,ind(k-1,2),ind(k-1,3));
                    exp_inv_temp = log(rand)/-lam;
                end
                cost_vec(obin,3) = cost(lam,min(s,net_size)+1);
                spell = exp_inv_temp;
                cum_spell = spell + cum_spell;
            end
        end
    lag = deathind(n)+1;
    end
end
    
