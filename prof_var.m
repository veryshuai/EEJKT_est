function res = prof_var(full_sim_dat)
% This function creates a smoothed 3D plot of profit variance conditional
% productivity phi and success probability theta

    % Read in data
    dat = full_sim_dat{10};
    prods = full(dat(:,1));
    mu_f = full(dat(:,2));
    
    % Get unique values
    prod_list = full(unique(prods));
    mu_f_list = full(unique(mu_f));
     
    % Get Productivity Coef of Variation
    sale_f_mat = full_sim_dat{4};
    rev_coef_var = zeros(size(prod_list,1),size(mu_f_list,1));
    for kind = 1:size(prod_list)
        for jind = 1:size(mu_f_list)
            rev_temp = sale_f_mat(:,prods == prod_list(kind) & mu_f == mu_f_list(jind)); 
            if max(rev_temp(:) > 0)
                rev_coef_var(kind,jind) = sqrt(var(rev_temp(:) > 0)) / mean(rev_temp(:) > 0);
            end
        end
    end
    
    % Plot
    surf(rev_coef_var);    

end %end function
