% This script snippet reads productivities into the event matrix, and success probabilities into a separate matrix, for a single firm

prod_init = find(rand<cum_erg_pp,1,'first');%initial productivity drawn out of ergodic distribution 
ind(obin,2) = prod_init(1,1); %read initial productivity into even matrix
[succ_prob(obin,1),succ_prob(obin,2),succ_prob(obin,3),succ_prob(obin,4)] = update_succ_probs(succ_params); %draw success probability

lag = 1;
time = 0;
while time < TT
    exp_inv_temp = log(rand)/-d;
    time = time + exp_inv_temp; %time of death
    [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
    % COMMENTED SECTION ALLOWS PRODUCTIVITY TO MOVE AROUND FOR SINGLE FIRM, CURRENTLY DISABLED
    % death = find([0,cum_sp_p(j,lag:end)]'<time,1,'last'); %index of death
    % scalar_temp = death(1,1)-1; 
    % if scalar_temp ~= 0 
    %     for k = lag:lag+scalar_temp-1
    %         ind(obin,:) = [cum_sp_p(j,k),...
    %             ind(obin-1,2)+1-2*(rand<.5*(1+(ind(obin-1,2)-Phi_size-1)/Phi_size)),...
    %             0,0,-1,-1,-1,-1,zeros(1,2),-1,-1];
    %         [obin,max_mat_violation,vio] = step(obin,mms,max_mat_violation); if vio == 1 breakflag = 1; break; end
    %     end
    %     lag = scalar_temp+1;
    % end
    ind(obin,:) = [time,find(rand<cum_erg_pp,1,'first'),0,0,0,0,0,0,zeros(1,2),-1,-1];
    [succ_prob(obin,1),succ_prob(obin,2),succ_prob(obin,3),succ_prob(obin,4)] = update_succ_probs(succ_params); %draw success probability
    deathmat(obin,1) = 1; %new firm
end
