%This snippet adds shipment dummies to the sh matrices

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
                    cost_vec(obin,4) = F;
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
                    cost_vec(obin,2) = F;
                    sh_val_f(obin,2) = ind(k,7); %copy identity into sh_val
                end
                tot_time = tot_time + next_ship_time;
            end
        end
    end
end
