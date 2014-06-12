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
    
    % Get foreign and domestic sales data
    sale_f_mat = full_sim_dat{4}; %foreign data
    sale_h_mat = full_sim_dat{2}; %domestic data
    
    % Get Productivity Coef of Variation
    rev_coef_var = zeros(size(prod_list,1),size(mu_f_list,1)); %preallocate
    for kind = 1:size(prod_list) %loop through productivities
        for jind = 1:size(mu_f_list) %loop through thetas
            rev_temp = sale_f_mat(:,prods == prod_list(kind) & mu_f == mu_f_list(jind)); %only the firms with exactly that theta
            dom_rev  = sale_h_mat(:,prods == prod_list(kind) & mu_f == mu_f_list(jind)); %domestic sales for the same firms
            display(size(rev_temp));
            if max(rev_temp(:) > 0) %coef of var doesn't exist if mean = 0
                 %only firms with either exports or domestic sales are counted 
                rev_coef_var(kind,jind) = sqrt(var(rev_temp(:) > 0 | dom_rev(:) > 0)) / mean(rev_temp(:) > 0 | dom_rev(:) > 0);
            end
        end
    end
    
    % Plot surface
    figure(1);
    surf(rev_coef_var);    
    set(gca,'XTick',mu_f_list)      % set Xticks
    set(gca,'YTick',prod_list)      % set Yticks
    
end %end function
