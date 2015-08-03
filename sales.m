function [] = sales(scale_f,scale_h,de,maxc,Z,Phi,X_h,X_f,cf_num,increase,TT,S_old,t)
%this function takes scale factors, elasticity of demand, and the simulated
%state transition vector and returns sales per period between events for
%each firm.

shockyear = 30; %year of the policy shock

Z_big = [-inf,Z'];

for j = 1:S_old

        % Load written files
        load(sprintf('/gpfs/home/dcj138/scratch/temp_data/temp_%d_%d.mat', j, t))
        
        sale_h = cell(size(st_ind_cont));
        sale_f = cell(size(st_ind_cont));
        
        for k = 1:size(st_ind_cont,1)
            
            sale_h{k} = exp(scale_h)*exp(Phi(st_ind_cont{k}(:,2))).^(de-1).*exp(X_h(st_ind_cont{k}(:,3))).*sum(exp(Z_big(ds{k}(:,1:maxc)+1)).*sh{k}(:,1:maxc),2); 

            scale_f_vec = exp(scale_f); % case=not counterfactual

            if cf_num == 3
                % Insert macro shock
                shock_time = find(st_ind_cont{k}(:,1) < shockyear,1,'last'); %find shock time
                scale_f_vec = ones(size(st_ind_cont{k}(:,1))) * exp(scale_f); % read in shocked scale_f
                if isempty(shock_time) == 0
                    scale_f_vec(1:shock_time) = exp(scale_f - log(increase)); % read in unshocked scale_f
                end
            end

            sale_f{k} = scale_f_vec.*exp(Phi(st_ind_cont{k}(:,2))).^(de-1).*exp(X_f(st_ind_cont{k}(:,4))).*sum(exp(Z_big(ds{k}(:,maxc+1:2*maxc)+1)).*sh{k}(:,maxc+1:2*maxc),2);  
        end

        % Old variable names
        sale_f_cont = sale_f;
        sale_h_cont = sale_h;

        % Save to disk
        save(sprintf('/gpfs/home/dcj138/scratch/temp_data/temp_%d_%d.mat', j, t),'st_ind_cont' ,'ds' ,'sh' ,'sh_val_h' ,'sh_val_f' ,'cost_vec' ,'succ_prob','t','S','sale_h_cont','sale_f_cont')


    end

end
